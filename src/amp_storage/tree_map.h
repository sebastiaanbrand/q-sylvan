#ifndef TREE_MAP_H
#define TREE_MAP_H

#include "flt.h"

#ifdef __cplusplus
extern "C" {
#endif


void *tree_map_create(unsigned long max_size, double tolerance);
void tree_map_free(void *map);
int tree_map_find_or_put(void *dbs, fl_t val, unsigned long *ret);
int tree_map_find_or_put2(void *dbs, complex_t *v, unsigned long *ret);
fl_t *tree_map_get(void *dbs, unsigned long ref);
complex_t tree_map_get2(void *dbs, unsigned long ref);
unsigned long tree_map_num_entries(void *dbs);
double tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // TREE_MAP
