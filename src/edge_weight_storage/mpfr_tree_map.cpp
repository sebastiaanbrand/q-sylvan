#include "mpfr_tree_map.h"

#include <map>

// float "equality" tolerance
static mpfr_t TOLERANCE; // default 1e-14
static double TOLERANCE_double;


static bool mpfr_close(mpfr_srcptr x, mpfr_srcptr y) 
{
    if (TOLERANCE_double == 0.0) {
        return mpfr_eq(x, y, MPFR_PREC);
    }
    else {
        mpfr_t diff;
        mpfr_init2(diff, MPFR_PREC);
        mpfr_sub(diff, x, y, DEFAULT_RND);
        int cmp = mpfr_cmpabs(TOLERANCE, diff); // 1 if |a| > |b|, 0 if |a| == |b|, -1 if |a| < |b|
        mpfr_clear(diff); // not sure if we need to manually clear these
        return (cmp >= 0);
    }
}

static bool mpfr_custom_comp(mpfr_srcptr x, mpfr_srcptr y)
{
    if (mpfr_close(x, y)) return 0;
    else return mpfr_less_p(x, y);
}


struct custom_comparator
{
    bool operator() ( __mpfr_struct x, __mpfr_struct y ) const {
    	return mpfr_custom_comp(&x, &y);
    }
};


// map <key, val>
typedef std::map<__mpfr_struct, unsigned int, custom_comparator> map_mpfr_to_int_t;
typedef std::map<unsigned int, __mpfr_struct> map_int_to_mpfr_t;


// (two way) map to assign unique ints to mpfr vals, using std::map (RB trees)
struct mpfr_tree_map_s {
    map_mpfr_to_int_t* mpfr_to_int;
    map_int_to_mpfr_t* int_to_mpfr;
    unsigned int max_size;
    unsigned int entries;
};


// create()
mpfr_tree_map_t *
mpfr_tree_map_create(unsigned int ms, double tol)
{
    mpfr_tree_map_t *map = (mpfr_tree_map_t *) calloc(1, sizeof(mpfr_tree_map_t));
    mpfr_init2(TOLERANCE, MPFR_PREC);
    mpfr_set_d(TOLERANCE, tol, DEFAULT_RND);
    TOLERANCE_double = tol;
    map->max_size = ms;
    map->entries = 0;
    map->mpfr_to_int = new map_mpfr_to_int_t;
    map->int_to_mpfr = new map_int_to_mpfr_t;
    return map;
}


// free()
void
mpfr_tree_map_free(mpfr_tree_map_t *map)
{
    delete map->mpfr_to_int;
    delete map->int_to_mpfr;
    free(map);
}


// find_or_put()
int
mpfr_tree_map_find_or_put(mpfr_tree_map_t *map, mpfr_t val, unsigned int *ret)
{
    map_mpfr_to_int_t *f2i = map->mpfr_to_int;
    map_int_to_mpfr_t *i2f = map->int_to_mpfr;

    // sigh... 
    __mpfr_struct val_cp;
    mpfr_init2(&val_cp, MPFR_PREC);
    mpfr_set(&val_cp, val, DEFAULT_RND);

    // look for double
    auto it = f2i->find(val_cp);
    if ( it == f2i->end() ) {
        // check if space
        if (map->entries > map->max_size-2) {
            printf("AMP map full\n");
            return -1;
        }

        // if key not found, insert new pair
        std::pair<__mpfr_struct, unsigned int> f2i_pair = std::make_pair(val_cp, map->entries);
        f2i->insert(f2i_pair);
        *ret = map->entries;

        // also insert in reverse table for lookup purposes
        std::pair<unsigned int, __mpfr_struct> i2f_pair = std::make_pair(map->entries, val_cp);
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
mpfr_ptr
mpfr_tree_map_get(mpfr_tree_map_t *map, unsigned int ref)
{
    return &(map->int_to_mpfr->find(ref)->second);   
}

/*
// size()
unsigned int
tree_map_size(tree_map_t *map)
{
    return map->entries;
}

// get_tolerance
double
tree_map_get_tolerance()
{
    return TOLERANCE;
}
*/


