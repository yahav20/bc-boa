#ifndef POOL_H
#define POOL_H

/* Forward-declare to avoid the circular include chain
 * pool.h -> node.h -> include.h -> heap.h (which needs node types). */
struct snode;

struct snode* new_node(void);
void          pool_reset(void);
void          pool_free_all(void);  /* release all allocated chunks (pool_bench only) */

#endif
