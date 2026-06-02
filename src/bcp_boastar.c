/////////////////////////////////////////////////////////////////////
// Bounded-Cost BOA* (BCP-BOA*)
// Faithful implementation of Skyler et al., SoCS 2022
/////////////////////////////////////////////////////////////////////

#include "heap.h"
#include "node.h"
#include "include.h"
#include "boastar.h"
#include "bcp_boastar.h"
#include "pool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>

extern gnode* graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start;
extern gnode* start_state;
extern gnode* goal_state;

extern unsigned long long int stat_expansions;
extern unsigned long long int stat_generated;
extern unsigned long long int stat_created;
extern unsigned long long int stat_recycled;
extern unsigned long long int minf_solution;
extern unsigned solutions[MAX_SOLUTIONS][2];
extern unsigned nsolutions;
extern unsigned stat_pruned;

snode* recycled_nodes[MAX_RECYCLE];
int next_recycled = 0;

/*
 * Global search-round counter.  Incremented once per bc_boastar() call.
 * Each gnode resets its dominance fields lazily on first access this round,
 * eliminating the O(num_gnodes) init loop that would otherwise touch
 * hundreds of MB per call.
 */
static unsigned g_bc_version = 0;

/*
 * Pareto-front storage for Min/Max/Avg dominance checks.
 * Stored OUTSIDE gnode to keep gnode ~64 bytes (cache-friendly).
 *
 * MAX_PARETO_PER_NODE=200 covers the paper's largest reported BAY POF of
 * 147 solutions.  Each row is state_id * MAX_PARETO_PER_NODE.
 * g_pareto_store / g_pareto_count are allocated lazily on the first
 * bc_boastar() call after the graph has been loaded.
 */
#define MAX_PARETO_PER_NODE 200

typedef struct { unsigned g1; unsigned g2; } ParetoPoint;

static ParetoPoint *g_pareto_store = NULL; /* [num_gnodes * MAX_PARETO_PER_NODE] */
static unsigned    *g_pareto_count = NULL; /* [num_gnodes] */

static void ensure_pareto_store(void) {
    if (g_pareto_store != NULL) return;
    g_pareto_store = (ParetoPoint*)calloc((size_t)num_gnodes * MAX_PARETO_PER_NODE,
                                           sizeof(ParetoPoint));
    g_pareto_count = (unsigned*)calloc(num_gnodes, sizeof(unsigned));
    if (!g_pareto_store || !g_pareto_count) {
        fprintf(stderr, "bc_boastar: failed to allocate Pareto store (%u nodes)\n",
                num_gnodes);
        exit(1);
    }
}

static inline void lazy_reset(gnode *n, unsigned state_id) {
    n->version  = g_bc_version;
    n->gmin     = LARGE;
    n->g1min    = LARGE;
    if (g_pareto_store)
        g_pareto_count[state_id] = 0;
}

/*
 * Returns 1 if path (g1,g2) to state_id is dominated by a previously
 * expanded path recorded for this round.
 *
 * Lex1 / Lex2 : O(1) — single scalar compare, identical to BOA*.
 * Min/Max/Avg : O(|front|) — scan the true non-dominated Pareto front.
 */
static int is_dominated(unsigned state_id, unsigned g1, unsigned g2, OrderingFunction ord) {
    gnode *n = &graph_node[state_id];

    if (n->version != g_bc_version) return 0; /* fresh this round */

    if (ord == ORDER_LEX1) return (g2 >= n->gmin);
    if (ord == ORDER_LEX2) return (g1 >= n->g1min);

    /* Min / Max / Avg */
    ParetoPoint *front = g_pareto_store + (size_t)state_id * MAX_PARETO_PER_NODE;
    unsigned cnt = g_pareto_count[state_id];
    for (unsigned i = 0; i < cnt; i++)
        if (front[i].g1 <= g1 && front[i].g2 <= g2) return 1;
    return 0;
}

/*
 * Record that (g1,g2) was expanded at state_id.
 * For Min/Max/Avg: maintains a true (non-dominated) Pareto front,
 * removing entries dominated by (g1,g2) before appending the new point.
 * MUST be called exactly once per expansion, after is_dominated returned 0.
 */
static void record_dominance(unsigned state_id, unsigned g1, unsigned g2, OrderingFunction ord) {
    gnode *n = &graph_node[state_id];

    if (n->version != g_bc_version)
        lazy_reset(n, state_id);

    if (ord == ORDER_LEX1) { if (g2 < n->gmin)  n->gmin  = g2; return; }
    if (ord == ORDER_LEX2) { if (g1 < n->g1min) n->g1min = g1; return; }

    /* Min / Max / Avg: maintain true Pareto front. */
    ParetoPoint *front = g_pareto_store + (size_t)state_id * MAX_PARETO_PER_NODE;
    unsigned *cnt = &g_pareto_count[state_id];

    unsigned write = 0;
    for (unsigned i = 0; i < *cnt; i++) {
        if (front[i].g1 >= g1 && front[i].g2 >= g2) continue; /* dominated */
        front[write++] = front[i];
    }
    *cnt = write;
    if (*cnt < MAX_PARETO_PER_NODE) {
        front[*cnt].g1 = g1;
        front[*cnt].g2 = g2;
        (*cnt)++;
    }
    /* If every slot holds a truly non-dominated point and we're still full,
     * we omit the new point.  The search stays correct — just slightly less
     * aggressive at pruning for states with >200 non-dominated paths. */
}

/* Timeout support: set bc_timeout_ms > 0 before calling bc_boastar().
 * bc_timed_out is set to 1 if the search aborted due to timeout. */
double bc_timeout_ms = 0.0;
int    bc_timed_out  = 0;

/* -----------------------------------------------------------------------
 * Compute the key for a node given an ordering function.
 * SEL_LEX must be resolved to LEX1 or LEX2 before calling.
 * ----------------------------------------------------------------------- */
static double compute_key(unsigned f1, unsigned f2, OrderingFunction ord,
                           unsigned long long min1, unsigned long long max1,
                           unsigned long long min2, unsigned long long max2)
{
    double nf1, nf2, fmin, fmax, favg;

    switch (ord) {
    case ORDER_LEX1:
        return (double)f1 * (double)BASE + (double)f2;

    case ORDER_LEX2:
        return (double)f2 * (double)BASE + (double)f1;

    case ORDER_MIN:
        nf1  = (double)(f1 - min1) / (double)(max1 - min1);
        nf2  = (double)(f2 - min2) / (double)(max2 - min2);
        fmin = nf1 < nf2 ? nf1 : nf2;
        fmax = nf1 > nf2 ? nf1 : nf2;
        return fmin * (double)BASE + fmax;

    case ORDER_MAX:
        nf1  = (double)(f1 - min1) / (double)(max1 - min1);
        nf2  = (double)(f2 - min2) / (double)(max2 - min2);
        fmin = nf1 < nf2 ? nf1 : nf2;
        fmax = nf1 > nf2 ? nf1 : nf2;
        return fmax * (double)BASE + fmin;

    case ORDER_AVG:
        nf1  = (double)(f1 - min1) / (double)(max1 - min1);
        nf2  = (double)(f2 - min2) / (double)(max2 - min2);
        fmin = nf1 < nf2 ? nf1 : nf2;
        favg = (nf1 + nf2) / 2.0;
        return favg * (double)BASE + fmin;

    default:
        return (double)f1 * (double)BASE + (double)f2;
    }
}

/* -----------------------------------------------------------------------
 * Compute map extremes: min1, min2, max2, max1.
 * Also sets h1 and h2 for ALL graph nodes via backward Dijkstra.
 * ----------------------------------------------------------------------- */
void compute_map_extremes(unsigned long long *min1, unsigned long long *max2,
                           unsigned long long *max1, unsigned long long *min2)
{
    if (backward_dijkstra(1)) *min1 = start_state->h1;
    if (backward_dijkstra(2)) *min2 = start_state->h2;

    /* Forward Lex(c1,c2): key = (c1 << 32) | c2 */
    for (unsigned i = 0; i < num_gnodes; ++i)
        graph_node[i].key = 0xFFFFFFFFFFFFFFFFULL;
    emptyheap_dij();
    start_state->key = 0;
    insertheap_dij(start_state);

    while (topheap_dij() != NULL) {
        gnode* n = popheap_dij();
        if (n == goal_state) break;
        unsigned long long c1 = n->key >> 32;
        unsigned long long c2 = n->key & 0xFFFFFFFFULL;
        for (short d = 1; d < (short)(adjacent_table[n->id][0] * 3); d += 3) {
            gnode* succ = &graph_node[adjacent_table[n->id][d]];
            unsigned long long new_key =
                ((c1 + adjacent_table[n->id][d + 1]) << 32) |
                 (c2 + adjacent_table[n->id][d + 2]);
            if (succ->key > new_key) {
                succ->key = new_key;
                insertheap_dij(succ);
            }
        }
    }
    *max2 = goal_state->key & 0xFFFFFFFFULL;

    /* Forward Lex(c2,c1): key = (c2 << 32) | c1 */
    for (unsigned i = 0; i < num_gnodes; ++i)
        graph_node[i].key = 0xFFFFFFFFFFFFFFFFULL;
    emptyheap_dij();
    start_state->key = 0;
    insertheap_dij(start_state);

    while (topheap_dij() != NULL) {
        gnode* n = popheap_dij();
        if (n == goal_state) break;
        unsigned long long c2 = n->key >> 32;
        unsigned long long c1 = n->key & 0xFFFFFFFFULL;
        for (short d = 1; d < (short)(adjacent_table[n->id][0] * 3); d += 3) {
            gnode* succ = &graph_node[adjacent_table[n->id][d]];
            unsigned long long new_key =
                ((c2 + adjacent_table[n->id][d + 2]) << 32) |
                 (c1 + adjacent_table[n->id][d + 1]);
            if (succ->key > new_key) {
                succ->key = new_key;
                insertheap_dij(succ);
            }
        }
    }
    *max1 = goal_state->key & 0xFFFFFFFFULL;
}

/* -----------------------------------------------------------------------
 * BCP-BOA* main search.
 * b1, b2       : absolute cost bounds.
 * ord          : ordering function (SEL_LEX must be pre-resolved).
 * min1..max2   : normalization extremes (for MIN/MAX/AVG orderings).
 * Returns 1 if a solution within bounds was found, 0 otherwise.
 * ----------------------------------------------------------------------- */
int bc_boastar(unsigned b1, unsigned b2, OrderingFunction ord,
               unsigned long long min1, unsigned long long max1,
               unsigned long long min2, unsigned long long max2)
{
    struct timeval bc_t0;
    gettimeofday(&bc_t0, NULL);

    next_recycled   = 0;
    nsolutions      = 0;
    stat_pruned     = 0;
    stat_expansions = 0;
    stat_generated  = 0;
    stat_recycled   = 0;
    stat_created    = 0;
    bc_timed_out    = 0;

    emptyheap();

    /* Allocate the external Pareto store once (after graph is loaded). */
    ensure_pareto_store();

    /*
     * Advance the global round counter.  Every gnode resets itself lazily
     * (in record_dominance) the first time it is accessed this round.
     * This replaces a slow O(num_gnodes) memset that would touch hundreds
     * of MB for each of the 90+ bc_boastar calls in a benchmark query.
     */
    g_bc_version++;

    if (start_state->h1 > b1 || start_state->h2 > b2)
        return 0;

    snode* root = new_node();
    ++stat_created;
    root->state      = start_state->id;
    root->g1         = 0;
    root->g2         = 0;
    root->searchtree = NULL;
    root->key        = compute_key(start_state->h1, start_state->h2,
                                   ord, min1, max1, min2, max2);
    insertheap(root);

    while (topheap() != NULL) {
        snode* n = popheap();

        /* Timeout: fire every 16384 total pops (dominated + non-dominated).
         * Checking only on non-dominated expansions fails when most pops are
         * dominated (stat_expansions barely moves, timeout never fires). */
        if (bc_timeout_ms > 0.0 &&
            ((stat_expansions + (unsigned long long)stat_pruned) & 0x3FFF) == 0) {
            struct timeval now;
            gettimeofday(&now, NULL);
            double elapsed = 1000.0 * (now.tv_sec  - bc_t0.tv_sec)
                           + 0.001  * (now.tv_usec - bc_t0.tv_usec);
            if (elapsed >= bc_timeout_ms) {
                bc_timed_out = 1;
                emptyheap();
                return 0;
            }
        }

        /* Dominance check at expansion */
        if (is_dominated(n->state, n->g1, n->g2, ord)) {
            stat_pruned++;
            if (next_recycled < MAX_RECYCLE)
                recycled_nodes[next_recycled++] = n;
            continue;
        }

        /* Not dominated — record permanently so future nodes can be pruned */
        record_dominance(n->state, n->g1, n->g2, ord);

        /* Goal test */
        if (n->state == (int)goal_state->id) {
            solutions[nsolutions][0] = n->g1;
            solutions[nsolutions][1] = n->g2;
            nsolutions++;
            return 1;
        }

        ++stat_expansions;

        for (short d = 1; d < (short)(adjacent_table[n->state][0] * 3); d += 3) {
            unsigned nsucc = adjacent_table[n->state][d];
            unsigned cost1 = adjacent_table[n->state][d + 1];
            unsigned cost2 = adjacent_table[n->state][d + 2];

            unsigned newg1 = n->g1 + cost1;
            unsigned newg2 = n->g2 + cost2;
            unsigned h1    = graph_node[nsucc].h1;
            unsigned h2    = graph_node[nsucc].h2;
            unsigned f1    = newg1 + h1;
            unsigned f2    = newg2 + h2;

            /* Bound pruning */
            if (f1 > b1 || f2 > b2) {
                stat_pruned++;
                continue;
            }

            /* Early dominance check:
             * - Lex1/Lex2: O(1) scalar compare on gmin/g1min — always do it.
             * - Min/Max/Avg: O(|front|) pareto scan; skip at generation time
             *   to avoid billions of scans on hard queries.  Dominated nodes
             *   will be caught at expansion time instead. */
            if ((ord == ORDER_LEX1 || ord == ORDER_LEX2) &&
                is_dominated(nsucc, newg1, newg2, ord)) {
                stat_pruned++;
                continue;
            }

            /* Heap overflow guard: insertheap() has no bounds check.
             * Drop this successor rather than writing past heap[]. */
            if (heapsize >= 39000000UL) {
                stat_pruned++;
                continue;
            }

            snode* succ;
            if (next_recycled > 0) {
                succ = recycled_nodes[--next_recycled];
                stat_recycled++;
            } else {
                succ = new_node();
                ++stat_created;
            }

            succ->state      = nsucc;
            succ->g1         = newg1;
            succ->g2         = newg2;
            succ->searchtree = n;
            succ->key        = compute_key(f1, f2, ord, min1, max1, min2, max2);
            stat_generated++;
            insertheap(succ);
        }
    }

    return 0;
}

/* -----------------------------------------------------------------------
 * Select the 5 paper pivots from the full POF.
 * pof[i][0] = c1 cost, pof[i][1] = c2 cost of POF solution i.
 * Pivots are returned as normalised (nc1, nc2) ∈ [0,1]².
 * FTL = (0,1), FBR = (1,0) (the Lex extreme solutions).
 * MD, TL, BR are selected from the actual POF.
 * ----------------------------------------------------------------------- */
void select_pivots(unsigned (*pof)[2], unsigned npof,
                   unsigned long long min1, unsigned long long max1,
                   unsigned long long min2, unsigned long long max2,
                   double pivots[5][2])
{
    double range1 = (double)(max1 - min1);
    double range2 = (double)(max2 - min2);

    pivots[PIVOT_FTL][0] = 0.0;  pivots[PIVOT_FTL][1] = 1.0;
    pivots[PIVOT_FBR][0] = 1.0;  pivots[PIVOT_FBR][1] = 0.0;

    if (npof == 0 || range1 == 0.0 || range2 == 0.0) {
        pivots[PIVOT_MD][0] = 0.5; pivots[PIVOT_MD][1] = 0.5;
        pivots[PIVOT_TL][0] = 0.25; pivots[PIVOT_TL][1] = 0.75;
        pivots[PIVOT_BR][0] = 0.75; pivots[PIVOT_BR][1] = 0.25;
        return;
    }

    /* MD: POF solution minimising |nc1 - nc2| */
    int md_idx = 0;
    double md_diff = 1e18;
    for (unsigned i = 0; i < npof; ++i) {
        double nc1 = (pof[i][0] - min1) / range1;
        double nc2 = (pof[i][1] - min2) / range2;
        double diff = fabs(nc1 - nc2);
        if (diff < md_diff) { md_diff = diff; md_idx = (int)i; }
    }
    double md_nc1 = (pof[md_idx][0] - min1) / range1;
    double md_nc2 = (pof[md_idx][1] - min2) / range2;
    pivots[PIVOT_MD][0] = md_nc1;
    pivots[PIVOT_MD][1] = md_nc2;

    /* TL: closest POF solution to midpoint(FTL, MD) */
    double tl_mid0 = md_nc1 / 2.0;
    double tl_mid1 = (1.0 + md_nc2) / 2.0;
    int tl_idx = 0;
    double tl_dist = 1e18;
    for (unsigned i = 0; i < npof; ++i) {
        double nc1 = (pof[i][0] - min1) / range1;
        double nc2 = (pof[i][1] - min2) / range2;
        double d = (nc1-tl_mid0)*(nc1-tl_mid0) + (nc2-tl_mid1)*(nc2-tl_mid1);
        if (d < tl_dist) { tl_dist = d; tl_idx = (int)i; }
    }
    pivots[PIVOT_TL][0] = (pof[tl_idx][0] - min1) / range1;
    pivots[PIVOT_TL][1] = (pof[tl_idx][1] - min2) / range2;

    /* BR: closest POF solution to midpoint(FBR, MD) */
    double br_mid0 = (1.0 + md_nc1) / 2.0;
    double br_mid1 = md_nc2 / 2.0;
    int br_idx = 0;
    double br_dist = 1e18;
    for (unsigned i = 0; i < npof; ++i) {
        double nc1 = (pof[i][0] - min1) / range1;
        double nc2 = (pof[i][1] - min2) / range2;
        double d = (nc1-br_mid0)*(nc1-br_mid0) + (nc2-br_mid1)*(nc2-br_mid1);
        if (d < br_dist) { br_dist = d; br_idx = (int)i; }
    }
    pivots[PIVOT_BR][0] = (pof[br_idx][0] - min1) / range1;
    pivots[PIVOT_BR][1] = (pof[br_idx][1] - min2) / range2;
}

static const char* ordering_name(OrderingFunction ord) {
    switch (ord) {
    case ORDER_LEX1:    return "Lex1";
    case ORDER_LEX2:    return "Lex2";
    case ORDER_SEL_LEX: return "Sel-Lex";
    case ORDER_MIN:     return "Min";
    case ORDER_MAX:     return "Max";
    case ORDER_AVG:     return "Average";
    default:            return "?";
    }
}

static const char* pivot_name(PivotType p) {
    switch (p) {
    case PIVOT_FTL: return "FTL";
    case PIVOT_TL:  return "TL";
    case PIVOT_MD:  return "MD";
    case PIVOT_BR:  return "BR";
    case PIVOT_FBR: return "FBR";
    default:        return "?";
    }
}

/* -----------------------------------------------------------------------
 * Single-query wrapper (used by main_bc_boastar.c).
 * ----------------------------------------------------------------------- */
void call_bc_boastar(OrderingFunction ord, PivotType pivot_type, int zone)
{
    struct timeval t_start, t_setup_end, t_search_end;
    unsigned long long min1 = LARGE, min2 = LARGE;
    unsigned long long max1 = LARGE, max2 = LARGE;

    initialize_parameters();
    gettimeofday(&t_start, NULL);

    compute_map_extremes(&min1, &max2, &max1, &min2);

    printf("Map Extremes -> FTL: (%llu, %llu)  FBR: (%llu, %llu)\n",
           min1, max2, max1, min2);

    if (max1 == min1 || max2 == min2) {
        printf("Trivial: single-cost-dimension POF — using extreme solution.\n");
        printf("Solution: (%llu, %llu)\n", min1, min2);
        return;
    }

    boastar();

    unsigned npof = nsolutions;
    unsigned (*pof)[2] = malloc(npof * sizeof(*pof));
    if (!pof) { fprintf(stderr, "malloc failed for POF\n"); exit(1); }
    memcpy(pof, solutions, npof * sizeof(*pof));

    printf("BOA* found %u POF solutions.\n", npof);

    double pivots[5][2];
    select_pivots(pof, npof, min1, max1, min2, max2, pivots);
    free(pof);

    double px = pivots[pivot_type][0];
    double py = pivots[pivot_type][1];

    double b_bar1, b_bar2;
    switch (zone) {
    case 2: b_bar1=(3.0*px+1.0)/4.0; b_bar2=(3.0*py+1.0)/4.0; break;
    case 3: b_bar1=(px+1.0)/2.0;     b_bar2=(py+1.0)/2.0;     break;
    case 4: default:
            b_bar1=(px+3.0)/4.0;     b_bar2=(py+3.0)/4.0;     break;
    }

    unsigned b1 = (unsigned)(b_bar1*(double)(max1-min1)+(double)min1);
    unsigned b2 = (unsigned)(b_bar2*(double)(max2-min2)+(double)min2);

    OrderingFunction effective_ord = ord;
    if (ord == ORDER_SEL_LEX)
        effective_ord = (b_bar1 <= b_bar2) ? ORDER_LEX1 : ORDER_LEX2;

    printf("Pivot: %s (%.4f, %.4f)  Zone: %d  Bounds: b1=%u, b2=%u\n",
           pivot_name(pivot_type), px, py, zone, b1, b2);
    printf("Ordering: %s%s\n", ordering_name(ord),
           ord == ORDER_SEL_LEX
               ? (effective_ord == ORDER_LEX1 ? " -> Lex1" : " -> Lex2")
               : "");

    gettimeofday(&t_setup_end, NULL);

    bc_boastar(b1, b2, effective_ord, min1, max1, min2, max2);

    gettimeofday(&t_search_end, NULL);

    double setup_ms  = 1000.0*(t_setup_end.tv_sec  - t_start.tv_sec)
                     + 0.001 *(t_setup_end.tv_usec - t_start.tv_usec);
    double search_ms = 1000.0*(t_search_end.tv_sec  - t_setup_end.tv_sec)
                     + 0.001 *(t_search_end.tv_usec - t_setup_end.tv_usec);

    printf("\nBCP-BOA* Results\n");
    printf("----------------\n");
    printf("Start: %lld  Goal: %lld\n", start_state->id+1, goal_state->id+1);
    if (nsolutions > 0)
        printf("Solution: (%u, %u)\n", solutions[0][0], solutions[0][1]);
    else
        printf("No solution within bounds.\n");
    printf("Setup (heuristic+POF) time (ms): %.3f\n", setup_ms);
    printf("BCP-BOA* search time      (ms): %.3f\n", search_ms);
    printf("Total time                (ms): %.3f\n", setup_ms+search_ms);
    printf("States generated: %llu\n", stat_generated);
    printf("States expanded:  %llu\n", stat_expansions);
    printf("States pruned:    %u\n",   stat_pruned);
}
