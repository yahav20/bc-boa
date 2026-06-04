/*
 * benchmark.c — BCP-BOA* Benchmark
 * Replicates the experimental setup of Skyler et al., SoCS 2022.
 *
 * Usage:
 *   ./benchmark <map_file> <num_queries> <timeout_sec>
 *
 * Example:
 *   ./benchmark Maps/BAY-road-d.txt 50 30
 *
 * Query file is located automatically:
 *   "Maps/BAY-road-d.txt" -> "Queries/BAY-queries"
 *
 * For each query the benchmark runs:
 *   1. Heuristic computation (backward Dijkstra x2 + Lex Dijkstra x2)
 *   2. Full BOA* to obtain the Pareto-Optimal Frontier
 *   3. BCP-BOA* for every combination of:
 *        ordering  : Lex1, Lex2, Sel-Lex, Min, Max, Average
 *        pivot     : FTL, TL, MD, BR, FBR
 *        zone      : 2, 3, 4, 5
 *      (Zone 5 has b_bar=(1,1); all pivots give the same bounds so only
 *       one pivot is run for zone 5.)
 *
 * Output: per-query stats + summary tables mirroring paper Table 1.
 */

#include "include.h"
#include "boastar.h"
#include "bcp_boastar.h"
#include "graph.h"
#include "pool.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/* -----------------------------------------------------------------------
 * Externals needed from boastar.c / bcp_boastar.c
 * ----------------------------------------------------------------------- */
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
extern int next_recycled;          /* bcp_boastar.c recycle index */

/* -----------------------------------------------------------------------
 * Helpers
 * ----------------------------------------------------------------------- */
static inline double ms_diff(struct timeval a, struct timeval b) {
    return 1000.0*(b.tv_sec - a.tv_sec) + 0.001*(b.tv_usec - a.tv_usec);
}

/* -----------------------------------------------------------------------
 * Result accumulator for one (ordering, zone, pivot) combination.
 * ----------------------------------------------------------------------- */
typedef struct {
    double exp_sum;    /* sum of expansions across queries  */
    double gen_sum;    /* sum of generated                   */
    double prune_sum;  /* sum of pruned                      */
    double time_sum;   /* sum of search time (ms)            */
    double bbar1_sum;  /* sum of normalised b1               */
    double bbar2_sum;  /* sum of normalised b2               */
    int    solved;     /* queries where solution was found   */
    int    timed_out;  /* queries that hit the timeout       */
    int    n;          /* queries attempted                  */
    double error_sum;  /* sum of BS quality error (ORDER_BS only) */
    int    error_n;    /* queries where BS found a solution  */
} Acc;

/* Layout: [ordering 0..5][zone_idx 0..3][pivot_idx 0..4]
 * zone_idx: 0=Z2, 1=Z3, 2=Z4, 3=Z5   (Z5 uses only pivot_idx=0) */
#define N_ORD   7
#define N_ZONES 4
#define N_PIV   5

static Acc results[N_ORD][N_ZONES][N_PIV];

/* BOA* aggregate stats */
static double boa_exp_sum    = 0;
static double boa_gen_sum    = 0;
static double boa_sol_sum    = 0;
static double setup_time_sum = 0;   /* heuristic + BOA* time */
static int    n_queries_run  = 0;
static int    n_queries_skip = 0;

/* -----------------------------------------------------------------------
 * Zone formula: given normalised pivot (px,py) and zone, return b_bar.
 * ----------------------------------------------------------------------- */
static void zone_bounds(int zone, double px, double py,
                         double *bb1, double *bb2)
{
    switch (zone) {
    case 2: *bb1=(3.0*px+1.0)/4.0; *bb2=(3.0*py+1.0)/4.0; break;
    case 3: *bb1=(px+1.0)/2.0;     *bb2=(py+1.0)/2.0;     break;
    case 4: *bb1=(px+3.0)/4.0;     *bb2=(py+3.0)/4.0;     break;
    case 5: *bb1=1.0;               *bb2=1.0;               break;
    default:*bb1=(px+3.0)/4.0;     *bb2=(py+3.0)/4.0;
    }
}

/* -----------------------------------------------------------------------
 * Reset all node-pool and heap state between queries.
 * ----------------------------------------------------------------------- */
static void reset_for_query(void) {
    pool_free_all();     /* release all chunks so each query's peak is independent */
    boastar_pool_reset();
    next_recycled = 0;
}

/* -----------------------------------------------------------------------
 * Print one horizontal separator line of given width.
 * ----------------------------------------------------------------------- */
static void hline(int w) {
    for (int i=0;i<w;i++) putchar('-');
    putchar('\n');
}

/* -----------------------------------------------------------------------
 * Print the summary table for one metric.
 * metric[ord][zi][pi] = the per-combination average to display.
 * For Z5 only pi=0 is used (one column).
 * ----------------------------------------------------------------------- */
static void print_table(const char *title, double metric[N_ORD][N_ZONES][N_PIV],
                         Acc acc[N_ORD][N_ZONES][N_PIV],
                         const char *fmt)
{
    static const char *ord_names[] =
        {"Lex1","Lex2","Sel-Lex","Min","Max","Average","BS"};
    static const char *piv_names[] = {"FTL","TL","MD","BR","FBR"};
    static const int   zones[]     = {2,3,4,5};

    /* Column header indices: Z2×5 + Z3×5 + Z4×5 + Z5×1 = 16 columns */
    int ncols = N_PIV*3 + 1; /* 16 */
    int col_w = 8;
    int lbl_w = 10;

    printf("\n--- %s ---\n", title);

    /* Header row */
    printf("%-*s", lbl_w, "Ordering");
    for (int zi=0; zi<3; zi++)
        for (int pi=0; pi<N_PIV; pi++)
            printf(" Z%d:%-*s", zones[zi], col_w-4, piv_names[pi]);
    printf(" %-*s\n", col_w, "Z5:Any");

    /* b1 header row */
    printf("%-*s", lbl_w, "b1(avg)");
    for (int zi=0; zi<3; zi++)
        for (int pi=0; pi<N_PIV; pi++) {
            double avg = acc[0][zi][pi].n > 0
                         ? acc[0][zi][pi].bbar1_sum / acc[0][zi][pi].n : 0;
            printf(" %*.2f", col_w, avg);
        }
    {
        double avg = acc[0][3][0].n > 0
                     ? acc[0][3][0].bbar1_sum / acc[0][3][0].n : 1.0;
        printf(" %*.2f\n", col_w, avg);
    }

    /* b2 header row */
    printf("%-*s", lbl_w, "b2(avg)");
    for (int zi=0; zi<3; zi++)
        for (int pi=0; pi<N_PIV; pi++) {
            double avg = acc[0][zi][pi].n > 0
                         ? acc[0][zi][pi].bbar2_sum / acc[0][zi][pi].n : 0;
            printf(" %*.2f", col_w, avg);
        }
    {
        double avg = acc[0][3][0].n > 0
                     ? acc[0][3][0].bbar2_sum / acc[0][3][0].n : 1.0;
        printf(" %*.2f\n", col_w, avg);
    }

    hline(lbl_w + ncols*(col_w+1));

    /* Data rows */
    for (int oi=0; oi<N_ORD; oi++) {
        printf("%-*s", lbl_w, ord_names[oi]);
        for (int zi=0; zi<3; zi++) {
            for (int pi=0; pi<N_PIV; pi++) {
                if (acc[oi][zi][pi].n > 0)
                    printf(" %*.1f", col_w, metric[oi][zi][pi]);
                else
                    printf(" %*s", col_w, "N/A");
            }
        }
        /* Z5: only one column (pivot 0) */
        if (acc[oi][3][0].n > 0)
            printf(" %*.1f\n", col_w, metric[oi][3][0]);
        else
            printf(" %*s\n", col_w, "N/A");
    }
    (void)fmt;
}

/* -----------------------------------------------------------------------
 * main
 * ----------------------------------------------------------------------- */
int main(int argc, char **argv)
{
    if (argc < 4 || argc > 5) {
        fprintf(stderr,
            "Usage: %s <map_file> <num_queries> <timeout_sec> [start_query]\n"
            "  map_file     : e.g. Maps/BAY-road-d.txt\n"
            "  num_queries  : how many queries to run in this chunk\n"
            "  timeout_sec  : seconds allowed per ordering per query\n"
            "  start_query  : 0-indexed offset into query file (default 0)\n",
            argv[0]);
        return 1;
    }

    const char *map_file   = argv[1];
    int         num_qry    = atoi(argv[2]);
    double      timeout_s  = atof(argv[3]);
    double      timeout_ms = timeout_s * 1000.0;
    int         start_query = (argc >= 5) ? atoi(argv[4]) : 0;

    /* ---- Derive query file name from map file name ---- */
    /* "Maps/BAY-road-d.txt"  ->  prefix "BAY"
     * "Maps/FLA-road-d.txt"  ->  prefix "FLA"  etc. */
    const char *base = strrchr(map_file, '/');
    base = (base != NULL) ? base + 1 : map_file;
    char prefix[64] = {0};
    const char *road = strstr(base, "-road");
    const char *dash = strchr(base, '-');
    const char *cut  = (road != NULL) ? road : dash;
    if (cut != NULL) {
        int len = (int)(cut - base);
        if (len > 63) len = 63;
        strncpy(prefix, base, len);
    } else {
        strncpy(prefix, base, 63);
    }

    char query_file[256];
    snprintf(query_file, sizeof(query_file), "Queries/%s-queries", prefix);

    /* ---- Load graph ---- */
    read_adjacent_table(map_file);
    new_graph();

    /* ---- Read queries ---- */
    FILE *qf = fopen(query_file, "r");
    if (!qf) {
        fprintf(stderr, "Cannot open query file: %s\n", query_file);
        return 1;
    }

    unsigned queries[50][2];
    int nq = 0;
    while (nq < 50 &&
           fscanf(qf, "%u %u", &queries[nq][0], &queries[nq][1]) == 2)
        nq++;
    fclose(qf);

    if (nq == 0) { fprintf(stderr, "No queries read.\n"); return 1; }
    if (start_query < 0 || start_query >= nq) {
        fprintf(stderr, "Error: start_query=%d out of range (file has %d queries).\n",
                start_query, nq);
        return 1;
    }
    if (start_query + num_qry > nq) {
        int avail = nq - start_query;
        fprintf(stderr, "Warning: only %d queries available from offset %d (requested %d).\n",
                avail, start_query, num_qry);
        num_qry = avail;
    }

    /* ---- Header ---- */
    printf("=================================================================\n");
    printf("BCP-BOA* Benchmark\n");
    printf("Map   : %s (%u nodes)\n", map_file, num_gnodes);
    printf("Queries: %d (file indices %d-%d)  |  Timeout/ordering: %.1f s\n",
           num_qry, start_query+1, start_query+num_qry, timeout_s);
    printf("Query file: %s\n", query_file);
    printf("=================================================================\n\n");

    /* ---- Zero accumulators ---- */
    memset(results, 0, sizeof(results));

    /* ---- Main query loop ---- */
    for (int qi = start_query; qi < start_query + num_qry; qi++) {
        start = queries[qi][0] - 1;   /* file is 1-indexed */
        goal  = queries[qi][1] - 1;

        reset_for_query();

        int chunk_i = qi - start_query + 1;
        printf("Processing query %3d/%d  (file #%d: %u -> %u) ...\n",
               chunk_i, num_qry, qi+1, queries[qi][0], queries[qi][1]);
        fflush(stdout);

        /* Reset gmin for boastar() — g1min/pareto are handled lazily by bc_boastar */
        for (unsigned i = 0; i < num_gnodes; i++)
            graph_node[i].gmin = LARGE;

        initialize_parameters();

        /* --- Phase 1: heuristics + map extremes --- */
        struct timeval t0, t1;
        gettimeofday(&t0, NULL);

        unsigned long long min1 = LARGE, min2 = LARGE;
        unsigned long long max1 = LARGE, max2 = LARGE;
        compute_map_extremes(&min1, &max2, &max1, &min2);

        /* --- Phase 2: full BOA* (for POF and pivot selection) --- */
        stat_expansions = 0; stat_generated = 0;
        stat_pruned = 0; stat_created = 0; stat_recycled = 0;
        boastar();  /* also resets minf_solution internally */

        gettimeofday(&t1, NULL);
        double setup_ms = ms_diff(t0, t1);

        if (nsolutions == 0 || max1 == min1 || max2 == min2) {
            printf("Query %3d: skipped (trivial or no solution)\n", chunk_i);
            n_queries_skip++;
            continue;
        }
        
        printf("  -> BOA* baseline done: %u solutions. Running %d BCP combinations...\n",
               nsolutions, N_ORD * (3 * N_PIV + 1));
        fflush(stdout);

        double  boa_exp = (double)stat_expansions;
        double  boa_gen = (double)stat_generated;
        unsigned boa_sol = nsolutions;

        setup_time_sum += setup_ms;
        boa_exp_sum    += boa_exp;
        boa_gen_sum    += boa_gen;
        boa_sol_sum    += boa_sol;
        n_queries_run++;

        /* Save POF for pivot selection */
        unsigned npof = nsolutions;
        unsigned (*pof)[2] = malloc(npof * sizeof(*pof));
        if (!pof) { fprintf(stderr, "malloc failed\n"); return 1; }
        memcpy(pof, solutions, npof * sizeof(*pof));

        double pivots[5][2];
        select_pivots(pof, npof, min1, max1, min2, max2, pivots);
        /* pof is kept alive through the BCP phase for BS quality measurement */

        /* Free BOA* search nodes — they are no longer needed.
         * This lets the pool be reused for the BCP-BOA* phase,
         * keeping peak pool usage at max(BOA*, all-BCP) rather than
         * BOA* + all-BCP accumulated. */
        pool_reset();
        boastar_pool_reset();
        next_recycled = 0;

        /* --- Phase 3: BCP-BOA* for all combinations --- */
        /* Z2,Z3,Z4 × N_PIV pivots + Z5 × 1 pivot, for each of N_ORD orderings */
        const int total_combos = N_ORD * (3 * N_PIV + 1);
        int total_combos_run = 0;
        static const int zones[] = {2, 3, 4, 5};
        for (int zi = 0; zi < N_ZONES; zi++) {
            int zone   = zones[zi];
            int n_pivs = (zone == 5) ? 1 : N_PIV;

            for (int pi = 0; pi < n_pivs; pi++) {
                double px = pivots[pi][0];
                double py = pivots[pi][1];
                double bb1, bb2;
                zone_bounds(zone, px, py, &bb1, &bb2);
                unsigned b1 = (unsigned)(bb1*(double)(max1-min1)+(double)min1);
                unsigned b2 = (unsigned)(bb2*(double)(max2-min2)+(double)min2);

                for (int oi = 0; oi < N_ORD; oi++) {
                    OrderingFunction ord = (OrderingFunction)oi;
                    OrderingFunction eff = ord;
                    if (ord == ORDER_SEL_LEX)
                        eff = (bb1 <= bb2) ? ORDER_LEX1 : ORDER_LEX2;

                    bc_timeout_ms = timeout_ms;

                    pool_reset();
                    struct timeval ts, te;
                    gettimeofday(&ts, NULL);
                    bc_boastar(b1, b2, eff, min1, max1, min2, max2);
                    gettimeofday(&te, NULL);

                    double ms = ms_diff(ts, te);

                    Acc *a = &results[oi][zi][pi];
                    a->exp_sum   += (double)stat_expansions;
                    a->gen_sum   += (double)stat_generated;
                    a->prune_sum += (double)stat_pruned;
                    a->time_sum  += ms;
                    a->bbar1_sum += bb1;
                    a->bbar2_sum += bb2;
                    a->solved    += (nsolutions > 0 && !bc_timed_out);
                    a->timed_out += bc_timed_out;
                    a->n++;

                    /* BS quality: find nearest POF point by normalised Euclidean distance,
                     * then compute Error = (c1-p1)/p1 + (c2-p2)/p2 */
                    if (ord == ORDER_BS && nsolutions > 0 && !bc_timed_out && npof > 0) {
                        double range1 = (double)(max1 - min1);
                        double range2 = (double)(max2 - min2);
                        double c1 = (double)solutions[0][0];
                        double c2 = (double)solutions[0][1];
                        double best_dist = 1e18;
                        int    best_i    = 0;
                        for (unsigned pi2 = 0; pi2 < npof; pi2++) {
                            double d1 = (c1 - (double)pof[pi2][0]) / range1;
                            double d2 = (c2 - (double)pof[pi2][1]) / range2;
                            double dist = d1*d1 + d2*d2;
                            if (dist < best_dist) { best_dist = dist; best_i = (int)pi2; }
                        }
                        double p1  = (double)pof[best_i][0];
                        double p2  = (double)pof[best_i][1];
                        double err = (p1 > 0 ? (c1 - p1) / p1 : 0.0)
                                   + (p2 > 0 ? (c2 - p2) / p2 : 0.0);
                        a->error_sum += err;
                        a->error_n++;
                    }

                    total_combos_run++;
                    if (total_combos_run % 12 == 0 || total_combos_run == total_combos) {
                        printf("      ... %d/%d BCP runs completed\n", total_combos_run, total_combos);
                        fflush(stdout);
                    }
                }
            }
        }

        free(pof);   /* done with POF for this query */

        /* Per-query line */
        printf("Query %3d | setup %7.1f ms | BOA* exp=%llu gen=%llu sol=%u |",
               chunk_i, setup_ms,
               (unsigned long long)boa_exp,
               (unsigned long long)boa_gen, boa_sol);
        /* Show Lex1/MD/Z4 as a representative */
        {
            Acc *a = &results[ORDER_LEX1][2][PIVOT_MD];
            int  last = a->n - 1;
            (void)last;
            printf(" Lex1/MD/Z4: exp=%llu%s",
                   (unsigned long long)stat_expansions,
                   bc_timed_out ? " (TO)" : "");
        }
        printf("\n");
    }

    if (n_queries_run == 0) {
        printf("No valid queries completed.\n");
        return 0;
    }

    double N = (double)n_queries_run;

    /* ========================================================
     * Summary tables
     * ======================================================== */
    printf("\n=================================================================\n");
    printf("SUMMARY (%d/%d queries completed, %d skipped)\n",
           n_queries_run, num_qry, n_queries_skip);
    printf("=================================================================\n");

    printf("\n--- BOA* BASELINE (avg over %d queries) ---\n", n_queries_run);
    printf("  Setup (heuristic+BOA*) time : %7.2f ms\n", setup_time_sum/N);
    printf("  BOA* states expanded        : %7.0f\n",    boa_exp_sum/N);
    printf("  BOA* states generated       : %7.0f\n",    boa_gen_sum/N);
    printf("  POF solutions found         : %7.1f\n",    boa_sol_sum/N);

    /* Build average metric arrays for the three tables */
    double avg_exp  [N_ORD][N_ZONES][N_PIV];
    double avg_time [N_ORD][N_ZONES][N_PIV];
    double avg_prune[N_ORD][N_ZONES][N_PIV];
    double pct_ok   [N_ORD][N_ZONES][N_PIV];

    for (int oi=0;oi<N_ORD;oi++)
    for (int zi=0;zi<N_ZONES;zi++)
    for (int pi=0;pi<N_PIV;pi++) {
        Acc *a = &results[oi][zi][pi];
        if (a->n > 0) {
            avg_exp  [oi][zi][pi] = a->exp_sum   / a->n / 1000.0; /* thousands */
            avg_time [oi][zi][pi] = a->time_sum  / a->n;
            avg_prune[oi][zi][pi] = a->prune_sum / a->n / 1000.0;
            pct_ok   [oi][zi][pi] = 100.0 * a->solved / a->n;
        } else {
            avg_exp  [oi][zi][pi] = 0;
            avg_time [oi][zi][pi] = 0;
            avg_prune[oi][zi][pi] = 0;
            pct_ok   [oi][zi][pi] = 0;
        }
    }

    print_table("EXPANSIONS (avg, thousands)", avg_exp,  results, "%.1f");
    print_table("SEARCH TIME (avg, ms)",       avg_time, results, "%.2f");
    print_table("PRUNED (avg, thousands)",      avg_prune,results, "%.1f");
    print_table("SUCCESS RATE (%)",             pct_ok,  results, "%.0f");

    /* Generated nodes table */
    double avg_gen[N_ORD][N_ZONES][N_PIV];
    for (int oi=0;oi<N_ORD;oi++)
    for (int zi=0;zi<N_ZONES;zi++)
    for (int pi=0;pi<N_PIV;pi++) {
        Acc *a = &results[oi][zi][pi];
        avg_gen[oi][zi][pi] = a->n > 0 ? a->gen_sum / a->n / 1000.0 : 0;
    }
    print_table("GENERATED (avg, thousands)", avg_gen, results, "%.1f");

    /* ---- BS Solution Quality table ---- */
    static const char *piv_names_q[] = {"FTL","TL","MD","BR","FBR"};
    static const int   zones_q[]     = {2,3,4,5};
    printf("\n--- BS SOLUTION QUALITY: avg error = (c1-p1)/p1 + (c2-p2)/p2 ---\n");
    printf("  (p1,p2) = nearest POF point by normalised Euclidean distance)\n");
    printf("  0.00 = exact Pareto point found;  N/A = no solution\n\n");
    printf("%-8s", "Zone");
    for (int zi = 0; zi < 3; zi++)
        for (int pi = 0; pi < N_PIV; pi++)
            printf(" Z%d:%-5s", zones_q[zi], piv_names_q[pi]);
    printf(" %-8s\n", "Z5:Any");
    hline(8 + (3*N_PIV+1)*8);
    {
        int oi = (int)ORDER_BS;
        printf("%-8s", "BS");
        for (int zi = 0; zi < 3; zi++) {
            for (int pi = 0; pi < N_PIV; pi++) {
                Acc *a = &results[oi][zi][pi];
                if (a->error_n > 0)
                    printf(" %7.4f", a->error_sum / a->error_n);
                else
                    printf(" %7s", "N/A");
            }
        }
        {
            Acc *a = &results[oi][3][0];
            if (a->error_n > 0)
                printf(" %7.4f\n", a->error_sum / a->error_n);
            else
                printf(" %7s\n", "N/A");
        }
    }

    printf("\n");

    /* ---- Append raw sums to CSV for cross-chunk aggregation ---- */
    {
        char csv_file[256];
        snprintf(csv_file, sizeof(csv_file), "%s_chunks.csv", prefix);

        /* Write header only if file is new/empty */
        int write_header = 0;
        FILE *tf = fopen(csv_file, "r");
        if (!tf) { write_header = 1; } else { fclose(tf); }

        FILE *cf = fopen(csv_file, "a");
        if (!cf) {
            fprintf(stderr, "Warning: cannot write chunk CSV to %s\n", csv_file);
        } else {
            if (write_header) {
                fprintf(cf, "# boa_stats: type,map,start_q,n,boa_exp_sum,boa_gen_sum,boa_sol_sum,setup_time_sum\n");
                fprintf(cf, "# bcp_stats: type,map,start_q,ord_i,zone_i,piv_i,exp_sum,gen_sum,prune_sum,time_sum,bbar1_sum,bbar2_sum,solved,timed_out,n,error_sum,error_n\n");
            }
            fprintf(cf, "boa_stats,%s,%d,%d,%.0f,%.0f,%.0f,%.4f\n",
                    prefix, start_query, n_queries_run,
                    boa_exp_sum, boa_gen_sum, boa_sol_sum, setup_time_sum);
            for (int oi = 0; oi < N_ORD; oi++)
            for (int zi = 0; zi < N_ZONES; zi++)
            for (int pi = 0; pi < N_PIV; pi++) {
                Acc *a = &results[oi][zi][pi];
                if (a->n == 0) continue;
                fprintf(cf, "bcp_stats,%s,%d,%d,%d,%d,%.0f,%.0f,%.0f,%.6f,%.8f,%.8f,%d,%d,%d,%.10f,%d\n",
                        prefix, start_query, oi, zi, pi,
                        a->exp_sum, a->gen_sum, a->prune_sum, a->time_sum,
                        a->bbar1_sum, a->bbar2_sum,
                        a->solved, a->timed_out, a->n,
                        a->error_sum, a->error_n);
            }
            fclose(cf);
            printf("Chunk data appended to: %s\n", csv_file);
        }
    }

    return 0;
}
