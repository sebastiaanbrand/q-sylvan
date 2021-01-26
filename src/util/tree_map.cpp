#include "tree_map.h"

#include <map>
#include <math.h>

// TODO: store fl_t instead of doubles

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;

static inline bool double_close(double x, double y) 
{
    return (fabs(x - y) < TOLERANCE);
}

static bool double_compare(double x, double y)
{
    if (double_close(x, y)) return 0;
    else return (x < y);
}

struct custom_comparator
{
    bool operator() ( double x, double y ) const {
    	return double_compare(x,y);
    }
};

// map <key, val>
typedef std::map<double, unsigned int, custom_comparator> map_double_to_int_t;
typedef std::map<unsigned int, double> map_int_to_double_t;


// (two way) map to assign unique integers to doubles, using std::map (RB trees)
struct tree_map_s {
    map_double_to_int_t* double_to_int;
    map_int_to_double_t* int_to_double;
    unsigned int entries;
};


// tree_map_create()
tree_map_t *
tree_map_create(double tol)
{
    tree_map_t *map = (tree_map_t *) calloc(1, sizeof(tree_map_t));
    TOLERANCE = tol;
    map->entries = 0;
    map->double_to_int = new map_double_to_int_t;
    map->int_to_double = new map_int_to_double_t;
    return map;
}

// tree_map_free()
void
tree_map_free(tree_map_t *map)
{
    delete map->double_to_int;
    delete map->int_to_double;
    free(map);
}

// TODO: tree_map_find_or_put()
int
tree_map_find_or_put(tree_map_t *map, double val, unsigned int *ret)
{
    map_double_to_int_t *d2i = map->double_to_int;
    map_int_to_double_t *i2d = map->int_to_double;

    // look for double
    auto it = d2i->find(val);
    if ( it == d2i->end() ) {
        // if key not found, insert new pair
        std::pair<double, unsigned int> d2i_pair = std::make_pair(val, map->entries);
        d2i->insert(d2i_pair);
        *ret = map->entries;

        // also insert in reverse table for lookup purposes
        std::pair<unsigned int, double> i2d_pair = std::make_pair(map->entries, val);
        i2d->insert(i2d_pair);

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
double *
tree_map_get(tree_map_t *map, unsigned int ref)
{
    // both these syntactically beautiful statements should be equivalent
    //return &(*(map->int_to_double))[ref];
    return &(map->int_to_double->find(ref)->second);
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


