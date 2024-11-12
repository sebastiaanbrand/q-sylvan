#ifndef AMP_STORAGE_INTERFACE
#define AMP_STORAGE_INTERFACE

#include "flt.h"
#include "cmap.h"
#include "mpc_map.h"

typedef enum wgt_storage_backend {
    COMP_HASHMAP,
    n_wgt_storage_types,
    MPC_HASHMAP
} wgt_storage_backend_t;

// create(uint64_t size, double tolerance)
extern void * (*wgt_store_create)(uint64_t size, double tolerance);

// free
extern void (*wgt_store_free)(void *wgt_storage);

// find_or_put(void *dbs, void *v, int *ret)
extern int (*wgt_store_find_or_put)(const void *dbs, const void *v, uint64_t *ret);

// get(void *dbs, int ref)
extern void * (*wgt_store_get)(const void *dbs, const uint64_t ref);

// num entries(void *dbs)
extern uint64_t (*wgt_store_num_entries)(const void *dbs);

// get tolerance
extern double (*wgt_store_get_tol)();

void init_wgt_storage_functions(wgt_storage_backend_t backend);

#endif // AMP_STORAGE_INTERFACE
