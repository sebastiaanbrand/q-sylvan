#include "mpfr_tree_map.h"

#include <map>

// float "equality" tolerance

// TODO: mpfr
static long double TOLERANCE = 1e-14l;

/*
static inline bool flt_close(fl_t x, fl_t y) 
{
    return (flt_abs(x - y) < TOLERANCE);
}

static bool flt_compare(fl_t x, fl_t y)
{
    if (flt_close(x, y)) return 0;
    else return (x < y);
}

struct custom_comparator
{
    bool operator() ( fl_t x, fl_t y ) const {
    	return flt_compare(x,y);
    }
};
*/

// map <key, val>
typedef std::map<mpfr_t, unsigned int> map_mpfr_to_int_t; // custom_comparator
typedef std::map<unsigned int, mpfr_t> map_int_to_mpfr_t;


// (two way) map to assign unique integers to doubles, using std::map (RB trees)
struct mpfr_tree_map_s {
    map_mpfr_to_int_t* mpfr_to_int;
    map_int_to_mpfr_t* int_to_mpfr;
    unsigned int max_size;
    unsigned int entries;
};


// tree_map_create()
mpfr_tree_map_t *
mpfr_tree_map_create(unsigned int ms, double tol)
{
    mpfr_tree_map_t *map = (mpfr_tree_map_t *) calloc(1, sizeof(mpfr_tree_map_t));
    TOLERANCE = tol;
    map->max_size = ms;
    map->entries = 0;
    map->mpfr_to_int = new map_mpfr_to_int_t;
    map->int_to_mpfr = new map_int_to_mpfr_t;
    return map;
}


// tree_map_free()
void
mpfr_tree_map_free(mpfr_tree_map_t *map)
{
    delete map->mpfr_to_int;
    delete map->int_to_mpfr;
    free(map);
}

/*

// TODO: tree_map_find_or_put()
int
tree_map_find_or_put(tree_map_t *map, fl_t val, unsigned int *ret)
{
    map_flt_to_int_t *f2i = map->flt_to_int;
    map_int_to_flt_t *i2f = map->int_to_flt;

    // look for double
    auto it = f2i->find(val);
    if ( it == f2i->end() ) {
        // check if space
        if (map->entries > map->max_size-2) {
            printf("AMP map full\n");
            return -1;
        }

        // if key not found, insert new pair
        std::pair<fl_t, unsigned int> f2i_pair = std::make_pair(val, map->entries);
        f2i->insert(f2i_pair);
        *ret = map->entries;

        // also insert in reverse table for lookup purposes
        std::pair<unsigned int, fl_t> i2f_pair = std::make_pair(map->entries, val);
        i2f->insert(i2f_pair);

        map->entries++;
        return 0;
    }
    else {
        // if it is found, return the (unique) int corresponding to this double
        *ret = it->second;
        return 1;
    }
}

// tree_map_get()
fl_t *
tree_map_get(tree_map_t *map, unsigned int ref)
{
    // both these syntactically beautiful statements should be equivalent
    //return &(*(map->int_to_flt))[ref];
    return &(map->int_to_flt->find(ref)->second);
}

// tree_map_size()
unsigned int
tree_map_size(tree_map_t *map)
{
    return map->entries;
}

// tree_map_get_tolerance
double
tree_map_get_tolerance()
{
    return TOLERANCE;
}
*/


