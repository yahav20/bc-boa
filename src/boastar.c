/////////////////////////////////////////////////////////////////////
// Carlos Hernandez
// All rights reserved
/////////////////////////////////////////////////////////////////////

#include "heap.h"
#include "node.h"
#include "include.h"
#include "boastar.h"
#include "pool.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>

// Global variables for graph and search state

gnode* graph_node; // Array of graph nodes
unsigned num_gnodes; // Number of nodes in the graph
unsigned adjacent_table[MAXNODES][MAXNEIGH]; // Adjacency list for successors
unsigned pred_adjacent_table[MAXNODES][MAXNEIGH]; // Adjacency list for predecessors
unsigned goal, start; // Goal and start node indices

gnode* start_state; // Pointer to start node

gnode* goal_state; // Pointer to goal node
snode* start_node; // Pointer to the start search node

// Statistics and solution storage
unsigned long long int stat_expansions = 0; // Number of node expansions
unsigned long long int stat_generated = 0;  // Number of nodes generated
unsigned long long int minf_solution = LARGE; // Cost of best solution found
unsigned long long int stat_created = 0; // Number of snode structs created
unsigned long long int stat_recycled = 0; // recycled counter

unsigned solutions[MAX_SOLUTIONS][2]; // Stores found solutions (g1, g2)
unsigned nsolutions = 0; // Number of solutions found
unsigned stat_pruned = 0; // Number of pruned nodes

// Initializes start and goal pointers and resets percolation stats
void initialize_parameters() {
    start_state = &graph_node[start];
    goal_state = &graph_node[goal];
    stat_percolations = 0;
}

// Runs backward Dijkstra from the goal to compute heuristic values (h1 or h2)
// dim==1: cost heuristic, dim==2: time heuristic
int backward_dijkstra(int dim) {
    for (int i = 0; i < num_gnodes; ++i)
        graph_node[i].key = LARGE; // Reset all node keys
    emptyheap_dij(); // Clear Dijkstra heap
    goal_state->key = 0; // Goal has zero cost to itself
    insertheap_dij(goal_state); // Insert goal into heap

    while (topheap_dij() != NULL) {
        gnode* n;
        gnode* pred;
        short d;
        n = popheap_dij(); // Get node with smallest key
        if (dim == 1)
            n->h1 = n->key; // Set cost heuristic
        else
            n->h2 = n->key; // Set time heuristic
        ++stat_expansions;
        // For each predecessor of n
        for (d = 1; d < pred_adjacent_table[n->id][0] * 3; d += 3) {
            pred = &graph_node[pred_adjacent_table[n->id][d]];
            int new_weight = n->key + pred_adjacent_table[n->id][d + dim];
            if (pred->key > new_weight) {
                pred->key = new_weight;
                insertheap_dij(pred);
            }
        }
    }
    return 1;
}

/* File-static recycle pool — avoids the 8 MB stack array of the original code. */
static snode* boa_recycled_nodes[MAX_RECYCLE];
static int    boa_next_recycled = 0;

/* Called by the benchmark before each query to invalidate stale pool pointers. */
void boastar_pool_reset(void) {
    boa_next_recycled = 0;
}

// Main BOA* search algorithm
int boastar() {
    boa_next_recycled = 0;
    nsolutions = 0;
    minf_solution = LARGE;
    stat_pruned = 0;

    emptyheap(); // Clear open list

    // Initialize start node
    start_node = new_node();
    ++stat_created;
    start_node->state = start;
    start_node->g1 = 0; // Cost so far (dimension 1)
    start_node->g2 = 0; // Cost so far (dimension 2)
    start_node->key = 0; // Key for heap ordering
    start_node->searchtree = NULL; // No parent
    insertheap(start_node); // Add to open list

    stat_expansions = 0;

    // Main search loop
    while (topheap() != NULL) {
        snode* n = popheap(); // Get best node from open list
        short d;

        // Prune if dominated or not promising
        if (n->g2 >= graph_node[n->state].gmin || n->g2 + graph_node[n->state].h2 >= minf_solution) {
            stat_pruned++;
            if (boa_next_recycled < MAX_RECYCLE) {
                boa_recycled_nodes[boa_next_recycled++] = n;
            }
            continue;
        }

        graph_node[n->state].gmin = n->g2; // Update best g2 for this node

        // Goal test
        if (n->state == goal) {
            // Store solution
            solutions[nsolutions][0] = n->g1;
            solutions[nsolutions][1] = n->g2;
            nsolutions++;
            if (nsolutions >= MAX_SOLUTIONS) {
                printf("Maximum number of solutions reached, increase MAX_SOLUTIONS!\n");
                exit(1);
            }
            if (minf_solution > n->g2)
                minf_solution = n->g2; // Update best solution
            continue;
        }

        ++stat_expansions;

        // Expand successors
        for (d = 1; d < adjacent_table[n->state][0] * 3; d += 3) {
            snode* succ;
            double newk1, newk2, newkey;
            unsigned nsucc = adjacent_table[n->state][d]; // Successor node index
            unsigned cost1 = adjacent_table[n->state][d + 1]; // Cost for dim 1
            unsigned cost2 = adjacent_table[n->state][d + 2]; // Cost for dim 2

            unsigned newg1 = n->g1 + cost1;
            unsigned newg2 = n->g2 + cost2;
            unsigned h1 = graph_node[nsucc].h1;
            unsigned h2 = graph_node[nsucc].h2;

            // Prune if dominated or not promising
            if (newg2 >= graph_node[nsucc].gmin || newg2 + h2 >= minf_solution)
                continue;

            newk1 = newg1 + h1;
            newk2 = newg2 + h2;

            // Reuse pruned node if available
            if (boa_next_recycled > 0) {
                succ = boa_recycled_nodes[--boa_next_recycled];
                stat_recycled++;
            }
            else {
                succ = new_node();
                ++stat_created;
            }

            succ->state = nsucc;
            stat_generated++;

            newkey = newk1 * (double)BASE + newk2; // Key for heap ordering
            succ->searchtree = n; // Set parent
            succ->g1 = newg1;
            succ->g2 = newg2;
            succ->key = newkey;
            insertheap(succ); // Add to open list
        }
    }

    return nsolutions > 0; // Return true if any solution found
}

/* ------------------------------------------------------------------------------*/
// Wrapper to run BOA* and print results
void call_boastar() {
    float runtime, heuristic_runtime;
    struct timeval tstart, tend;
    unsigned long long min_cost;
    unsigned long long min_time;
    struct timeval compute_heuristic_time;
    initialize_parameters(); // Set up start/goal pointers

    gettimeofday(&tstart, NULL); // Start timer

    // Compute heuristics using backward Dijkstra
    if (backward_dijkstra(1))
        min_cost = start_state->h1;
    if (backward_dijkstra(2))
        min_time = start_state->h2;

//    compute_heuristic(HEURISTIC_PESSIMISTIC_DEGREE);

    gettimeofday(&compute_heuristic_time, NULL); // end heuristic timer
    // Run BOA*
    boastar();

    gettimeofday(&tend, NULL); // End timer
    runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    heuristic_runtime = 1.0 * (compute_heuristic_time.tv_sec - tstart.tv_sec) + 1.0 * (compute_heuristic_time.tv_usec - tstart.tv_usec) / 1000000.0;
    // Print results
    printf("BOA* Search Results\n");
    printf("-------------------\n");
    printf("Start Node: %lld\n", start_state->id + 1);
    printf("Goal Node: %lld\n", goal_state->id + 1);
    printf("Number of Solutions: %d\n", nsolutions);
    printf("Runtime (ms): %.3f\n", runtime * 1000);
    printf("Heuristic Runtime (ms): %.3f\n", heuristic_runtime * 1000);
    printf("BOA* Runtime (ms): %.3f\n", (runtime - heuristic_runtime) * 1000);
    printf("States Generated: %llu\n", stat_generated);
    printf("States Expanded: %llu\n", stat_expansions);
    printf("States Created: %llu\n", stat_created);
    printf("Recycled: %llu\n", stat_recycled);
    // Verify all solutions not dominated
    bool dominated = false;
    for (int i = 0; i < nsolutions; i++)
    {
        for (int j = 0 ; j < nsolutions; j++)
        {
            if (solutions[i][0] <= solutions[j][0] && solutions[i][1] <= solutions[j][1] && i != j)
            {
                dominated = true;
                // printf("Solution %d is dominated by solution %d\n", i, j);
                // printf("Solution %d: (%lld, %lld)\n", i, solutions[i][0], solutions[i][1]);
                // printf("Solution %d: (%lld, %lld)\n", j, solutions[j][0], solutions[j][1]);
                break;
            }
        }
    }
    printf("Dominated: %s\n", dominated ? "Yes" : "No");
}
