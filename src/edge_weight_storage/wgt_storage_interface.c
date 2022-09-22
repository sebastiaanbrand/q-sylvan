#include "wgt_storage_interface.h"

void init_wgt_storage_functions(wgt_storage_backend_t backend)
{

    switch (backend)
    {
    case COMP_HASHMAP:
        wgt_store_create      = &cmap_create;
        wgt_store_free        = &cmap_free;
        wgt_store_find_or_put = &cmap_find_or_put;
        wgt_store_get         = &cmap_get;
        wgt_store_num_entries = &cmap_count_entries;
        wgt_store_get_tol     = &cmap_get_tolerance;
        break;
    case REAL_HASHMAP:
        wgt_store_create      = &rmap_create;
        wgt_store_free        = &rmap_free;
        wgt_store_find_or_put = &rmap_find_or_put2;
        wgt_store_get         = &rmap_get2;
        wgt_store_num_entries = &rmap_count_entries;
        wgt_store_get_tol     = &rmap_get_tolerance;
        break;
    case REAL_TREE:
        wgt_store_create      = &tree_map_create;
        wgt_store_free        = &tree_map_free;
        wgt_store_find_or_put = &tree_map_find_or_put2;
        wgt_store_get         = &tree_map_get2;
        wgt_store_num_entries = &tree_map_num_entries;
        wgt_store_get_tol     = &tree_map_get_tolerance;
    default:
        break;
    }
}
