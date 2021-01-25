#include "tree_map.h"

#include <map>

// float "equality" tolerance
// TODO: maybe make the map_t a struct, with the std::map as one component, 
// and these static variables as other components. That allows multiple 
// instances if needed.
static long double TOLERANCE = 1e-14l;
static int entries = 0;

// map <key, val>
typedef std::map<double, unsigned int> double_map_t;


// TODO: change data type from void* to something clearer
void *
tree_map_create(double tolerance)
{
    TOLERANCE = tolerance;
    return (void *) (new double_map_t);
}

// tree_map_free()
void
tree_map_free(void *map)
{
    double_map_t *m = (double_map_t *)(map);
    delete m;
}

// tree_map_find_or_put()
int
tree_map_find_or_put(void *map, double val, unsigned int *ret)
{
    double_map_t *m = (double_map_t *)(map);

    // look for double
    auto it = m->find(val);
    if ( it == m->end() ) {
        // if key not found, insert new pair
        std::pair<unsigned int, double> new_pair = std::make_pair(val, entries);
        m->insert(new_pair);
        *ret = entries;
        entries++;
        return 0;
    }
    else {
        // if it is found, return the (unique) int corresponding to this double
        *ret = it->second;
        return 1;
    }

}

// TODO: tree_map_get()

// tree_map_size()
int
tree_map_size(void *map)
{
    double_map_t *m = (double_map_t *)(map);
    return m->size();
}

// tree_map_get_tolerance
double
tree_map_get_tolerance()
{
    return TOLERANCE;
}


