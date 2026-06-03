/* pool.c — malloc-based allocator used by the regular boa/bcp_boastar binaries.
 * Memory is not freed between runs (acceptable for single-query binaries).
 */
#include "pool.h"
#include "include.h"   /* full snode definition */
#include <stdlib.h>
#include <stdio.h>

snode* new_node(void) {
    snode* n = (snode*)malloc(sizeof(snode));
    if (!n) { fprintf(stderr, "new_node: malloc failed\n"); exit(1); }
    n->heapindex = 0;
    return n;
}

void pool_reset(void)    { /* no-op for malloc version */ }
void pool_free_all(void) { /* no-op for malloc version */ }
