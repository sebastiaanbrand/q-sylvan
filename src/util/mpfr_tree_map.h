#ifndef MPFR_TREE_MAP_H
#define MPFR_TREE_MAP_H

#include <gmp.h>
#include <mpfr.h>

#define MPFR_PREC 40
#define DEFAULT_RND MPFR_RNDN

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mpfr_tree_map_s mpfr_tree_map_t;

mpfr_tree_map_t *mpfr_tree_map_create(unsigned int max_size, double tolerance);
void mpfr_tree_map_free(mpfr_tree_map_t *map);
int mpfr_tree_map_find_or_put(mpfr_tree_map_t *map, mpfr_t *val, unsigned int *ret);
//fl_t *mpfr_tree_map_get(mpfr_tree_map_t *map, unsigned int ref);
//unsigned int mpfr_tree_map_size(mpfr_tree_map_t *map);
//double mpfr_tree_map_get_tolerance();



#ifdef __cplusplus
}
#endif


#endif // MPFR_TREE_MAP
