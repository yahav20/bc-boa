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

typedef enum {
    ORDER_LEX1    = 0,
    ORDER_LEX2    = 1,
    ORDER_SEL_LEX = 2,
    ORDER_MIN     = 3,
    ORDER_MAX     = 4,
    ORDER_AVG     = 5
} OrderingFunction;

typedef enum {
    PIVOT_FTL = 0,
    PIVOT_TL  = 1,
    PIVOT_MD  = 2,
    PIVOT_BR  = 3,
    PIVOT_FBR = 4
} PivotType;

/* Timeout support: set before calling bc_boastar(); check after. */
extern double bc_timeout_ms;
extern int    bc_timed_out;

/* Public API */
void call_bc_boastar(OrderingFunction ord, PivotType pivot, int zone);

void compute_map_extremes(unsigned long long *min1, unsigned long long *max2,
                           unsigned long long *max1, unsigned long long *min2);

int  bc_boastar(unsigned b1, unsigned b2, OrderingFunction ord,
                unsigned long long min1, unsigned long long max1,
                unsigned long long min2, unsigned long long max2);

void select_pivots(unsigned (*pof)[2], unsigned npof,
                   unsigned long long min1, unsigned long long max1,
                   unsigned long long min2, unsigned long long max2,
                   double pivots[5][2]);

#endif
