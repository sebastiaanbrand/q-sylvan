#include "mpreal_tree_map.h"

#include <map>

// float "equality" tolerance
static mpfr::mpreal TOLERANCE = 1e-14;


static bool complex_custom_comp(mpreal_complex x, mpreal_complex y)
{
    if ((mpfr::abs(x.r - y.r) < TOLERANCE)) {
        if ((mpfr::abs(x.i - y.i) < TOLERANCE)) {
            // both within tolerance
            return 0;
        }
        else {
            return (x.i < y.i);
        }
    }
    else {
        return (x.r < y.r);
    }
}

struct custom_comparator
{
    bool operator() (mpreal_complex x, mpreal_complex y) const {
    	return complex_custom_comp(x, y);
    }
};

// map <key, val>
typedef std::map<mpreal_complex, unsigned int, custom_comparator> map_mpreal_to_int_t; // custom_comparator
typedef std::map<unsigned int, mpreal_complex> map_int_to_mpreal_t;


// (two way) map to assign unique ints to mpfr vals, using std::map (RB trees)
struct mpreal_tree_map_s {
    map_mpreal_to_int_t* mpreal_to_int;
    map_int_to_mpreal_t* int_to_mpreal;
    unsigned int max_size;
    unsigned int entries;
};



// create()
mpreal_tree_map_t *
mpreal_tree_map_create(unsigned int ms, double tol)
{
    mpfr::mpreal::set_default_prec(MPREAL_PREC);
    mpreal_tree_map_t *map = (mpreal_tree_map_t *) calloc(1, sizeof(mpreal_tree_map_t));
    if (tol > 0) TOLERANCE = tol;
    map->max_size = ms;
    map->entries = 0;
    map->mpreal_to_int = new map_mpreal_to_int_t;
    map->int_to_mpreal = new map_int_to_mpreal_t;
    return map;
}


// free()
void
mpreal_tree_map_free(mpreal_tree_map_t *map)
{
    delete map->mpreal_to_int;
    delete map->int_to_mpreal;
    free(map);
}


// find_or_put()
int
mpreal_tree_map_find_or_put(mpreal_tree_map_t *map, mpreal_complex val, unsigned int *ret)
{
    map_mpreal_to_int_t *f2i = map->mpreal_to_int;
    map_int_to_mpreal_t *i2f = map->int_to_mpreal;

    
    // look for double
    auto it = f2i->find(val);
    if ( it == f2i->end() ) {
        // check if space
        if (map->entries > map->max_size-2) {
            printf("MPREAL AMP map full\n");
            return -1;
        }

        // if key not found, insert new pair
        std::pair<mpreal_complex, unsigned int> f2i_pair = std::make_pair(val, map->entries);
        f2i->insert(f2i_pair);
        *ret = map->entries;

        // also insert in reverse table for lookup purposes
        std::pair<unsigned int, mpreal_complex> i2f_pair = std::make_pair(map->entries, val);
        i2f->insert(i2f_pair);

        map->entries++;
        return 0;
    }
    else {
        // if it is found, return the (unique) int corresponding to this double
        *ret = it->second;
        return 1;
    }
   return 0;
}


// get()
mpreal_complex
mpreal_tree_map_get(mpreal_tree_map_t *map, unsigned int ref)
{
    return map->int_to_mpreal->find(ref)->second;   
}


// size()
unsigned int
mpreal_tree_map_size(mpreal_tree_map_t *map)
{
    return map->entries;
}

// get_tolerance
double
mpreal_tree_map_get_tolerance()
{
    return TOLERANCE.toDouble();
}


