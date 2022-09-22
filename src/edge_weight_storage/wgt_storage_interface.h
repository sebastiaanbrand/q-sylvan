#ifndef AMP_STORAGE_INTERFACE
#define AMP_STORAGE_INTERFACE

#include "flt.h"
#include "cmap.h"
#include "rmap.h"
#include "tree_map.h"

typedef enum wgt_storage_backend {
    COMP_HASHMAP, 
    REAL_HASHMAP, 
    REAL_TREE,
    n_backends
} wgt_storage_backend_t;

// create(uint64_t size, double tolerance)
void * (*wgt_store_create)();

// free
void (*wgt_store_free)();

// find_or_put(void *dbs, complex_t *v, int *ret)
int (*wgt_store_find_or_put)();

// get(void *dbs, int ref)
complex_t (*wgt_store_get)();

// num entries(void *dbs)
unsigned long (*wgt_store_num_entries)();

// get tolerance
double (*wgt_store_get_tol)();

void init_wgt_storage_functions(wgt_storage_backend_t backend);

#endif // AMP_STORAGE_INTERFACE
