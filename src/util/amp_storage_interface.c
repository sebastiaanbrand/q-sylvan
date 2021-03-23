#include <stdio.h>

#include "amp_storage_interface.h"
#include "cmap.h"
#include "rmap.h"
#include "tree_map.h" 

void init_amp_storage_functions()
{
    // create
    amp_store_create[COMP_HASHMAP] = &cmap_create;
    amp_store_create[REAL_HASHMAP] = &rmap_create;
    amp_store_create[REAL_TREE] = &tree_map_create;
    
    // free
    amp_store_free[COMP_HASHMAP] = &cmap_free;
    amp_store_free[REAL_HASHMAP] = &rmap_free;
    amp_store_free[REAL_TREE] = &tree_map_free;

    // get tolerance
    amp_store_get_tol[COMP_HASHMAP] = &cmap_get_tolerance;
    amp_store_get_tol[REAL_HASHMAP] = &rmap_get_tolerance;
    amp_store_get_tol[REAL_TREE] = &tree_map_get_tolerance;
}