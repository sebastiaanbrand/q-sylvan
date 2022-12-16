#include "tree_map.h"

#include <map>
#include <math.h>
#include <cstdio>

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;

static inline bool flt_close(fl_t x, fl_t y) 
{
    if (TOLERANCE == 0.0) {
        return x == y;
    }
    else {
        return (flt_abs(x - y) < TOLERANCE);
    }
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

// map <key, val>
typedef std::map<fl_t, unsigned long, custom_comparator> map_flt_to_int_t;
typedef std::map<unsigned long, fl_t> map_int_to_flt_t;


// (two way) map to assign unique integers to doubles, using std::map (RB trees)
typedef struct tree_map_s tree_map_t;
struct tree_map_s {
    map_flt_to_int_t* flt_to_int;
    map_int_to_flt_t* int_to_flt;
    uint64_t max_size;
    uint64_t entries;
};

static unsigned long
tree_map_pack_indices(const tree_map_t * map, unsigned long r, unsigned long i)
{
    int index_size = (int) ceil(log2(map->max_size));
    unsigned long res = ((r << index_size) | i);
    return res;
}

static void
tree_map_unpack_indices(const tree_map_t * rmap, unsigned long bundle, unsigned long *index_r, unsigned long *index_i)
{
    unsigned long index_size = (long) ceil(log2(rmap->max_size));
    *index_r = (bundle >> index_size);
    *index_i = (bundle & ((1<<index_size)-1));
}


// tree_map_create()
void *
tree_map_create(uint64_t ms, double tol)
{
    tree_map_t *map = (tree_map_t *) calloc(1, sizeof(tree_map_t));
    TOLERANCE = tol;
    map->max_size = ms;
    map->entries = 0;
    map->flt_to_int = new map_flt_to_int_t;
    map->int_to_flt = new map_int_to_flt_t;
    return (void *) map;
}

// tree_map_free()
void
tree_map_free(void *m)
{
    tree_map_t * map = (tree_map_t *) m;
    delete map->flt_to_int;
    delete map->int_to_flt;
    free(map);
}

// tree_map_find_or_put()
int
tree_map_find_or_put(const void *dbs, const fl_t val, uint64_t *ret)
{
    tree_map_t * map = (tree_map_t *) dbs;
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
        std::pair<fl_t, unsigned long> f2i_pair = std::make_pair(val, map->entries);
        f2i->insert(f2i_pair);
        *ret = map->entries;

        // also insert in reverse table for lookup purposes
        std::pair<unsigned long, fl_t> i2f_pair = std::make_pair(map->entries, val);
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

int
tree_map_find_or_put2(const void *dbs, const complex_t *v, uint64_t *ret)
{
    uint64_t index_r = 0, index_i = 0;
    int found_r = tree_map_find_or_put(dbs, v->r, &index_r);
    int found_i = tree_map_find_or_put(dbs, v->i, &index_i);
    if (found_r == -1 || found_i == -1) return -1;
    tree_map_t * map = (tree_map_t *) dbs;
    *ret = tree_map_pack_indices(map, index_r, index_i);
    return (found_r && found_i); // if at leat one not found, return 0
}

// tree_map_get()
fl_t *
tree_map_get(const void *dbs, const uint64_t ref)
{
    tree_map_t * map = (tree_map_t *) dbs;
    return &(map->int_to_flt->find(ref)->second); //return &(*(map->int_to_flt))[ref];
}

complex_t
tree_map_get2(const void *dbs, const uint64_t ref)
{
    complex_t res;
    unsigned long index_r, index_i;
    fl_t *r, *i;
    tree_map_t * map = (tree_map_t *) dbs;
    tree_map_unpack_indices(map, ref, &index_r, &index_i);
    r = tree_map_get(dbs, index_r);
    i = tree_map_get(dbs, index_i);
    res.r = *r;
    res.i = *i;
    return res;
}

// tree_map_num_entries()
uint64_t
tree_map_num_entries(const void *dbs)
{
    tree_map_t * map = (tree_map_t *) dbs;
    return map->entries;
}

// tree_map_get_tolerance
double
tree_map_get_tolerance()
{
    return TOLERANCE;
}


