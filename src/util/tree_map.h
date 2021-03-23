#ifndef TREE_MAP_H
#define TREE_MAP_H

#include "flt.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tree_map_s tree_map_t;

tree_map_t *tree_map_create(unsigned int max_size, double tolerance);
void tree_map_free(void *map);
int tree_map_find_or_put(tree_map_t *map, fl_t val, unsigned int *ret);
fl_t *tree_map_get(tree_map_t *map, unsigned int ref);
unsigned int tree_map_size(tree_map_t *map);
double tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // TREE_MAP
