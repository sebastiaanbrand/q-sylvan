#include "tree_map.h"

#include <map>

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;

typedef std::map<uint64_t, double> double_map_t;

extern "C" {

    // tree_map_create()
    void *
    tree_map_create(double tolerance)
    {
        TOLERANCE = tolerance;
        return (void *) (new double_map_t);
    }

    // tree_map_free()
    void
    tree_map_free(void * map)
    {
        double_map_t* m = (double_map_t*)(map);
        delete m;
    }

    // TODO: tree_map_find_or_put()
    // TODO: tree_map_get()
    // TODO: tree_map_count_entries()
    // TODO: tree_map_get_tolerance

}
