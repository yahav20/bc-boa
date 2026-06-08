/*
 * bcp_boastar.h — Public interface for BCP-BOA* (Bounded-Cost Bi-Objective A*)
 *
 * Includes the Budget Score (BS) ordering and relaxed scalar dominance,
 * novel contributions described in the accompanying research report.
 */

#ifndef BCP_BOASTARH
#define BCP_BOASTARH

#define MAX_SOLUTIONS 1000000
#define MAX_RECYCLE   1000000

#include "node.h"

extern gnode *graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start;

/*
 * Ordering functions that determine the open-list key and the dominance strategy.
 *
 * Lex1/Lex2/Sel-Lex use a single scalar for dominance — O(1) per check.
 * Min/Max/Average maintain a full Pareto front per node — O(|POF|) per check.
 * BS uses relaxed scalar dominance (minimum S value seen) — O(1) per check.
 */
typedef enum {
    ORDER_LEX1    = 0,  /* Primary: f1, tiebreak: f2.  Dominance: min g2 seen. */
    ORDER_LEX2    = 1,  /* Primary: f2, tiebreak: f1.  Dominance: min g1 seen. */
    ORDER_SEL_LEX = 2,  /* Selects Lex1 if b_bar1 <= b_bar2, else Lex2.
                         * Must be resolved before calling bc_boastar(). */
    ORDER_MIN     = 3,  /* Primary: min(nf1,nf2), tiebreak: max.  Pareto dominance. */
    ORDER_MAX     = 4,  /* Primary: max(nf1,nf2), tiebreak: min.  Pareto dominance. */
    ORDER_AVG     = 5,  /* Primary: avg(nf1,nf2), tiebreak: min.  Pareto dominance. */
    ORDER_BS      = 6   /* Budget Score: S = (F1²+F2²)/(F1+F2), F1=f1/b1, F2=f2/b2.
                         * Scalar dominance: min S value seen per node. */
} OrderingFunction;

/*
 * The five pivot points used to select the budget centre on the POF.
 * All coordinates are normalised to [0,1]² by the range of optimal costs.
 *
 *   FTL (0,1) — Far Top-Left:     Lex-optimal in dimension 1.
 *   TL        — Top-Left:         Closest POF point to midpoint(FTL, MD).
 *   MD        — Middle-Diagonal:  Closest POF point to the nc1=nc2 diagonal.
 *   BR        — Bottom-Right:     Closest POF point to midpoint(FBR, MD).
 *   FBR (1,0) — Far Bottom-Right: Lex-optimal in dimension 2.
 */
typedef enum {
    PIVOT_FTL = 0,
    PIVOT_TL  = 1,
    PIVOT_MD  = 2,
    PIVOT_BR  = 3,
    PIVOT_FBR = 4
} PivotType;

/* Timeout support: set bc_timeout_ms > 0 before calling bc_boastar().
 * bc_timed_out is set to 1 if the search aborted early. */
extern double bc_timeout_ms;
extern int    bc_timed_out;

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

/* Single-query CLI wrapper: runs heuristics, BOA* (for POF), then BCP-BOA*. */
void call_bc_boastar(OrderingFunction ord, PivotType pivot, int zone);

/* Compute the four Lex-extreme costs and set h1/h2 for all graph nodes.
 * Must be called once per query before bc_boastar(). */
void compute_map_extremes(unsigned long long *min1, unsigned long long *max2,
                           unsigned long long *max1, unsigned long long *min2);

/* BCP-BOA* search.  Returns 1 if a solution within (b1,b2) was found, else 0.
 * ord must not be ORDER_SEL_LEX; resolve it first. */
int  bc_boastar(unsigned b1, unsigned b2, OrderingFunction ord,
                unsigned long long min1, unsigned long long max1,
                unsigned long long min2, unsigned long long max2);

/* Derive the five pivot points from the full POF (from BOA*).
 * Fills pivots[PIVOT_*][0..1] with normalised (nc1, nc2) coordinates. */
void select_pivots(unsigned (*pof)[2], unsigned npof,
                   unsigned long long min1, unsigned long long max1,
                   unsigned long long min2, unsigned long long max2,
                   double pivots[5][2]);

#endif
