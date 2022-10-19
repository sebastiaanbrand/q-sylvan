#ifndef TREE_MAP_H
#define TREE_MAP_H

#include "flt.h"

#ifdef __cplusplus
extern "C" {
#endif


void *tree_map_create(uint64_t max_size, double tolerance);
void tree_map_free(void *map);
int tree_map_find_or_put(const void *dbs, const fl_t val, uint64_t *ret);
int tree_map_find_or_put2(const void *dbs, const complex_t *v, uint64_t *ret);
fl_t *tree_map_get(const void *dbs, const uint64_t ref);
complex_t tree_map_get2(const void *dbs, const uint64_t ref);
uint64_t tree_map_num_entries(const void *dbs);
double tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // TREE_MAP
