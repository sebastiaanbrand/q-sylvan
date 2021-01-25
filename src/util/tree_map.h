#ifndef TREE_MAP_H
#define TREE_MAP_H

#ifdef __cplusplus
extern "C" {
#endif

void *tree_map_create(double tolerance);
void tree_map_free(void *map);
int tree_map_find_or_put(void *map, double val, unsigned int *ret);
int tree_map_size(void *map);
double tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // TREE_MAP
