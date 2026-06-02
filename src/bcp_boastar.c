/////////////////////////////////////////////////////////////////////
// Bounded-Cost BOA* (BCP-BOA*) - Lex1 Ordering
/////////////////////////////////////////////////////////////////////

#include "heap.h"
#include "node.h"
#include "include.h"
#include "boastar.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <stdbool.h>

// Global variables from your original setup
extern gnode* graph_node; 
extern unsigned num_gnodes; 
extern unsigned adjacent_table[MAXNODES][MAXNEIGH]; 
extern unsigned pred_adjacent_table[MAXNODES][MAXNEIGH]; 
extern unsigned goal, start; 
extern gnode* start_state; 
extern gnode* goal_state; 

// (Make sure these are declared in your globals or boastar.h)
extern unsigned long long int stat_expansions; 
extern unsigned long long int stat_generated;  
extern unsigned long long int minf_solution; 
extern unsigned long long int stat_created; 
extern unsigned long long int stat_recycled; 
extern unsigned solutions[MAX_SOLUTIONS][2]; 
extern unsigned nsolutions; 
extern unsigned stat_pruned; 

// Global variables to avoid stack overflow (MAX_RECYCLE is large)
snode* recycled_nodes[MAX_RECYCLE]; 
int next_recycled = 0;

snode* bc_start_node;


/**
 * Computes the extreme points (ideal and nadir) of the Pareto Front.
 * Uses backward Dijkstra for individual minima and lexicographical forward 
 * searches to find the corresponding maximum values for each dimension.
 */
void compute_map_extremes(unsigned long long *min1, unsigned long long *max2, unsigned long long *max1, unsigned long long *min2) {
    // 1. Compute backward heuristics (individual minima)
    if (backward_dijkstra(1)) *min1 = start_state->h1;
    if (backward_dijkstra(2)) *min2 = start_state->h2;
    
    // 2. Dynamic calculation of max2 using forward Lex(c1, c2) search
    for (int i = 0; i < num_gnodes; ++i) graph_node[i].key = 0xFFFFFFFFFFFFFFFFULL;
    emptyheap_dij();
    start_state->key = 0;
    insertheap_dij(start_state);
    
    while (topheap_dij() != NULL) {
        gnode* n = popheap_dij();
        if (n == goal_state) break;
        
        // Extract costs using bitmasking (32-bit each)
        unsigned long long current_c1 = n->key >> 32;
        unsigned long long current_c2 = n->key & 0xFFFFFFFFULL;
        
        for (short d = 1; d < adjacent_table[n->id][0] * 3; d += 3) {
            gnode* succ = &graph_node[adjacent_table[n->id][d]];
            unsigned cost1 = adjacent_table[n->id][d + 1];
            unsigned cost2 = adjacent_table[n->id][d + 2];
            
            unsigned long long new_key = ((current_c1 + cost1) << 32) | (current_c2 + cost2);
            if (succ->key > new_key) {
                succ->key = new_key;
                insertheap_dij(succ);
            }
        }
    }
    *max2 = goal_state->key & 0xFFFFFFFFULL;

    // 3. Dynamic calculation of max1 using forward Lex(c2, c1) search
    for (int i = 0; i < num_gnodes; ++i) graph_node[i].key = 0xFFFFFFFFFFFFFFFFULL;
    emptyheap_dij();
    start_state->key = 0;
    insertheap_dij(start_state);
    
    while (topheap_dij() != NULL) {
        gnode* n = popheap_dij();
        if (n == goal_state) break;
        
        unsigned long long current_c2 = n->key >> 32;
        unsigned long long current_c1 = n->key & 0xFFFFFFFFULL;
        
        for (short d = 1; d < adjacent_table[n->id][0] * 3; d += 3) {
            gnode* succ = &graph_node[adjacent_table[n->id][d]];
            unsigned cost1 = adjacent_table[n->id][d + 1];
            unsigned cost2 = adjacent_table[n->id][d + 2];
            
            unsigned long long new_key = ((current_c2 + cost2) << 32) | (current_c1 + cost1);
            if (succ->key > new_key) {
                succ->key = new_key;
                insertheap_dij(succ);
            }
        }
    }
    *max1 = goal_state->key & 0xFFFFFFFFULL;
}

// The main Bounded-Cost BOA* algorithm
// b1, b2 are the absolute cost bounds
int bc_boastar(unsigned b1, unsigned b2) {
    next_recycled = 0;
    nsolutions = 0;
    stat_pruned = 0;
    stat_expansions = 0;
    stat_generated = 0;
    stat_recycled = 0;

    emptyheap(); 
    // 1. Reset gmin for all nodes before the search to prevent dominance corruption.
    // This ensures that previous search runs do not interfere with current pruning logic.
    for (unsigned i = 0; i < num_gnodes; ++i) {
        graph_node[i].gmin = LARGE; 
    }

    // Initial bound pruning for start state - checked before allocation to avoid leak
    if (start_state->h1 > b1 || start_state->h2 > b2) {
        return 0; 
    }

    bc_start_node = new_node();
    ++stat_created;
    bc_start_node->state = start;

    // 2. Fix Indexing Mismatch by using start_state->id
    bc_start_node->state = start_state->id;
    bc_start_node->g1 = 0; 
    bc_start_node->g2 = 0; 
    bc_start_node->key = 0; 

    // 3. Proper Lex1 Key Initialization using f-values
    unsigned f1_start = start_state->h1;
    unsigned f2_start = start_state->h2;
    bc_start_node->key = ((unsigned long long)f1_start << 32) | (unsigned long long)f2_start;
    bc_start_node->searchtree = NULL; 
    
    insertheap(bc_start_node); 

    while (topheap() != NULL) {
        snode* n = popheap(); 
        short d;

        // Prune if dominated (Standard BOA* check)
        if (n->g2 >= graph_node[n->state].gmin) {
            stat_pruned++;
            if (next_recycled < MAX_RECYCLE) {
                recycled_nodes[next_recycled++] = n; 
            }
            continue;
        }

        graph_node[n->state].gmin = n->g2; 

        // Goal Test - Consistent indexing check for the goal state
        if (n->state == goal_state->id) {
            solutions[nsolutions][0] = n->g1;
            solutions[nsolutions][1] = n->g2;
            nsolutions++;
            return 1; // Terminate immediately upon finding a valid path
        }

        ++stat_expansions;

        // Expand successors
        for (d = 1; d < adjacent_table[n->state][0] * 3; d += 3) {
            snode* succ;
            unsigned long long newkey;
            unsigned nsucc = adjacent_table[n->state][d]; 
            unsigned cost1 = adjacent_table[n->state][d + 1]; 
            unsigned cost2 = adjacent_table[n->state][d + 2]; 

            unsigned newg1 = n->g1 + cost1;
            unsigned newg2 = n->g2 + cost2;
            unsigned h1 = graph_node[nsucc].h1;
            unsigned h2 = graph_node[nsucc].h2;
            
            unsigned f1 = newg1 + h1;
            unsigned f2 = newg2 + h2;

            // Bounded-Cost Pruning: discard if f-values exceed budgets 
            if (f1 > b1 || f2 > b2) {
                stat_pruned++;
                continue;
            }

            // Prune if dominated
            if (newg2 >= graph_node[nsucc].gmin)
                continue;

            if (next_recycled > 0) {
                succ = recycled_nodes[--next_recycled];
                stat_recycled++;
            } else {
                succ = new_node();
                ++stat_created;
            }

            succ->state = nsucc;
            stat_generated++;

            // Use bit-shifting for key calculation to maintain precision and Lex1 order
            newkey = ((unsigned long long)f1 << 32) | (unsigned long long)f2; 
            succ->searchtree = n; 
            succ->g1 = newg1;
            succ->g2 = newg2;
            succ->key = newkey;
            insertheap(succ); 
        }
    }

    return nsolutions > 0; 
}

/* ------------------------------------------------------------------------------*/
// Wrapper to setup bounds, run BCP-BOA* and print results
/* ------------------------------------------------------------------------------*/
void call_bc_boastar() {
    float runtime, heuristic_runtime;
    struct timeval tstart, tend, compute_heuristic_time;
    
    unsigned long long min1 = LARGE, min2 = LARGE;
    unsigned long long max1 = LARGE, max2 = LARGE;
    unsigned b1, b2; // Final absolute bounds

    initialize_parameters(); 
    gettimeofday(&tstart, NULL); 

    // 1. Compute Pareto Front extremes
    compute_map_extremes(&min1, &max2, &max1, &min2);

    gettimeofday(&compute_heuristic_time, NULL); 
    printf("Map Extremes Found -> Point 1: (%llu, %llu) | Point 2: (%llu, %llu)\n", min1, max2, max1, min2);

    /**
     * 2. Define bounds based on the BCP strategy.
     * Example: Zone 4 with a pivot point in the middle of the Nadir/Ideal range.
     * Formula: b_bar = ((pivot + 3) / 4)
     */
    double pivot_x = 0.5; 
    double pivot_y = 0.5;
    
    // Zone 4 formula from the paper: b_bar = ((x+3)/4, (y+3)/4)
    double b_bar_1 = (pivot_x + 3.0) / 4.0;
    double b_bar_2 = (pivot_y + 3.0) / 4.0;
    
    // De-normalization to get absolute budget values
    b1 = (unsigned)(b_bar_1 * (max1 - min1) + min1);
    b2 = (unsigned)(b_bar_2 * (max2 - min2) + min2);

    printf("Executing BCP-BOA* (Lex1) with Dynamic Bounds: b1=%u, b2=%u\n", b1, b2);

    // 3. Run the search
    bc_boastar(b1, b2);

    gettimeofday(&tend, NULL); 
    
    runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    heuristic_runtime = 1.0 * (compute_heuristic_time.tv_sec - tstart.tv_sec) + 1.0 * (compute_heuristic_time.tv_usec - tstart.tv_usec) / 1000000.0;
    
    printf("\nBCP-BOA* Search Results\n");
    printf("-------------------\n");
    printf("Start Node: %lld\n", start_state->id + 1);
    printf("Goal Node: %lld\n", goal_state->id + 1);
    
    if (nsolutions > 0) {
        printf("Solution Found within bounds: (%u, %u)\n", solutions[0][0], solutions[0][1]);
    } else {
        printf("No solution found within the specified bounds.\n");
    }
    
    printf("Total Runtime (ms): %.3f\n", runtime * 1000);
    printf("Pre-computation Runtime (ms): %.3f\n", heuristic_runtime * 1000);
    printf("Pure Search Runtime (ms): %.3f\n", (runtime - heuristic_runtime) * 1000);
    printf("States Generated: %llu\n", stat_generated);
    printf("States Expanded: %llu\n", stat_expansions);
    printf("States Pruned: %u\n", stat_pruned);
}