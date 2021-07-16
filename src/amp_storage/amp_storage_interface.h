#ifndef AMP_STORAGE_INTERFACE
#define AMP_STORAGE_INTERFACE

#include "flt.h"
#include "cmap.h"
#include "rmap.h"
#include "tree_map.h"

typedef enum amp_storage_backend {
    COMP_HASHMAP, 
    REAL_HASHMAP, 
    REAL_TREE,
    n_backends
} amp_storage_backend_t;

// create(uint64_t size, double tolerance)
void * (*amp_store_create[n_backends])();

// free
void (*amp_store_free[n_backends])();

// find_or_put(void *dbs, complex_t *v, int *ret)
int (*amp_store_find_or_put[n_backends])();

// get(void *dbs, int ref)
complex_t (*amp_store_get[n_backends])();

// num entries(void *dbs)
unsigned long (*amp_store_num_entries[n_backends])();

// get tolerance
double (*amp_store_get_tol[n_backends])();

void init_amp_storage_functions();

#endif // AMP_STORAGE_INTERFACE
