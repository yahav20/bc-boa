#ifndef BOASTARH
#define BOASTARH

#define MAX_SOLUTIONS 1000000
#define MAX_RECYCLE   1000000

extern gnode *graph_node;
extern unsigned num_gnodes;
extern unsigned adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
extern unsigned goal, start; 

void call_boastar();
void initialize_parameters();
int  backward_dijkstra(int dim);
int  boastar();
void boastar_pool_reset(void);

#endif
