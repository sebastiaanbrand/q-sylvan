#ifndef MPREAL_TREE_MAP_H
#define MPREAL_TREE_MAP_H

#include <gmp.h>
#include <mpreal.h>

#define MPREAL_PREC 500 // number of bits (for mantissa I guess)

typedef struct mpreal_tree_map_s mpreal_tree_map_t;

typedef struct
{
   mpfr::mpreal r;
   mpfr::mpreal i;
} mpreal_complex;

mpreal_tree_map_t *mpreal_tree_map_create(unsigned int max_size, double tolerance);
void mpreal_tree_map_free(mpreal_tree_map_t *map);
int mpreal_tree_map_find_or_put(mpreal_tree_map_t *map, mpreal_complex val, unsigned int *ret);
mpreal_complex mpreal_tree_map_get(mpreal_tree_map_t *map, unsigned int ref);
unsigned int mpreal_tree_map_size(mpreal_tree_map_t *map);
double mpreal_tree_map_get_tolerance();



#endif // MPREAL_TREE_MAP_H
