/* pool_bench.c — dynamic-chunk arena allocator for the benchmark binary.
 *
 * Grows in CHUNK_NODES-sized slabs on demand.  pool_reset() resets only
 * the bump pointers so the same physical pages are reused for the next query.
 * Peak footprint = largest single-query allocation (BOA* or BCP phases).
 */
#include "pool.h"
#include "include.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef CHUNK_NODES
#define CHUNK_NODES 1000000   /* nodes per slab (~56 MB each) */
#endif
#ifndef MAX_CHUNKS
#define MAX_CHUNKS 256        /* allocated on demand; actual peak depends on query difficulty */
#endif

static snode  *chunks[MAX_CHUNKS];
static int     n_chunks  = 0;
static int     cur_chunk = 0;
static unsigned cur_idx  = 0;

static void alloc_chunk(void) {
    if (n_chunks >= MAX_CHUNKS) {
        fprintf(stderr,
            "pool_bench: exceeded %d chunks (%d M nodes). "
            "Recompile with -DMAX_CHUNKS=<larger>.\n",
            MAX_CHUNKS, MAX_CHUNKS);
        exit(1);
    }
    chunks[n_chunks] = (snode*)malloc((size_t)CHUNK_NODES * sizeof(snode));
    if (!chunks[n_chunks]) {
        fprintf(stderr, "pool_bench: malloc failed for chunk %d\n", n_chunks);
        exit(1);
    }
    n_chunks++;
}

snode* new_node(void) {
    /* First call: no chunks yet */
    if (n_chunks == 0) {
        alloc_chunk();
        cur_chunk = 0;
        cur_idx   = 0;
    } else if (cur_idx >= (unsigned)CHUNK_NODES) {
        /* Current chunk full: advance to next, allocating if necessary */
        cur_chunk++;
        cur_idx = 0;
        if (cur_chunk >= n_chunks)
            alloc_chunk();
    }

    snode *n = &chunks[cur_chunk][cur_idx++];
    n->heapindex = 0;
    return n;
}

void pool_reset(void) {
    /* Rewind to start of first chunk; keep all slabs allocated for reuse */
    cur_chunk = 0;
    cur_idx   = 0;
}

void pool_free_all(void) {
    /* Free every allocated chunk and reset all counters.
     * Call between queries so each query's peak is independent. */
    for (int i = 0; i < n_chunks; i++) {
        free(chunks[i]);
        chunks[i] = NULL;
    }
    n_chunks  = 0;
    cur_chunk = 0;
    cur_idx   = 0;
}
