#ifndef TREE_MAP_H
#define TREE_MAP_H

#include "flt.h"

#ifdef __cplusplus
extern "C" {
#endif


void *tree_map_create(unsigned int max_size, double tolerance);
void tree_map_free(void *map);
int tree_map_find_or_put(void *dbs, fl_t val, unsigned int *ret);
int tree_map_find_or_put2(void *dbs, complex_t *v, unsigned int *ret);
fl_t *tree_map_get(void *dbs, unsigned int ref);
complex_t tree_map_get2(void *dbs, unsigned int ref);
unsigned int tree_map_size(void *dbs);
double tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // TREE_MAP
