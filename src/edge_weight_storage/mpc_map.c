#include "mpc_map.h"

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "atomics.h"
#include "fast_hash.h"
#include "util.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256

// how many "blocks" of 64 bits for a single table entry
#define entry_size (sizeof(mpc_ptr)/8)

typedef union {
    mpc_ptr         c;
    uint64_t        d[entry_size];
} bucket_t;

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;
static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

/**
\typedef Lockless hastable database.
*/
typedef struct mpc_map_s mpc_map_t;
struct mpc_map_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
    printf("\n");
}


int
mpc_map_find_or_put(const void *dbs, const void *_v, uint64_t *ret)
{
    mpc_ptr v = (mpc_ptr)_v;
    mpc_map_t *mpc_map = (mpc_map_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    // TODO: compute hash based on values (not on pointers)
    uint32_t hash = 0;
    //uint32_t hash  = SuperFastHash(&round_v, sizeof(complex_t), 0);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

    // Insert/lookup `v`
    for (unsigned int c = 0; c < mpc_map->threshold; c++) {
        uint64_t            ref = hash & mpc_map->mask;
        uint64_t            line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &mpc_map->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    // write backwards (overwrite bucket->d[0] last)
                    for (int k = entry_size-1; k >= 0; k--) {
                        atomic_write (&bucket->d[k], val->d[k]);
                    }
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            // 4. Bucket contains some complex value, check if close to `v`
            mpc_ptr in_table = (mpc_ptr) bucket;
            if (false) {
                // TODO: decide if values equal (or close enough)
                *ret = ref;
                return 1;
            }

            // If unsuccessful, try next
            ref += 1;
            ref = ref == line_end ? line_end - CACHE_LINE_SIZE : ref;
        }
        hash += prime << CACHE_LINE;
    }
    // amplitude table full, unable to add
    return -1;
}

void *
mpc_map_get(const void *dbs, const uint64_t ref)
{
    mpc_map_t *mpc_map = (mpc_map_t *) dbs;
    return &(mpc_map->table[ref].c);
}

uint64_t
mpc_map_count_entries(const void *dbs)
{
    mpc_map_t *mpc_map = (mpc_map_t *) dbs;
    uint64_t entries = 0;
    for (unsigned int c = 0; c < mpc_map->size; c++) {
        if (mpc_map->table[c].d[0] != EMPTY)
            entries++;
    }
    return entries;
}

void
mpc_map_print_bitvalues(const void *dbs, const uint64_t ref)
{
    mpc_map_t *mpc_map = (mpc_map_t *) dbs;
    bucket_t *b = (bucket_t *) mpc_map_get(mpc_map, ref);
    printf("%016" PRIu64, b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016" PRIu64, b->d[k]);
    }
}

void *
mpc_map_create(uint64_t size, double tolerance)
{
    mpc_map_t  *mpc_map = calloc (1, sizeof(mpc_map_t));
    mpc_map->size = size;
    mpc_map->mask = mpc_map->size - 1;
    mpc_map->table = calloc (mpc_map->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < mpc_map->size; c++) {
        mpc_map->table[c].d[0] = EMPTY;
    }
    mpc_map->threshold = mpc_map->size / 100;
    mpc_map->threshold = min(mpc_map->threshold, 1ULL << 16);
    mpc_map->seen_0 = 0;
    return (void *) mpc_map;
}

void
mpc_map_free(void *dbs)
{
    mpc_map_t * mpc_map = (mpc_map_t *) dbs;
    free (mpc_map->table);
    free (mpc_map);
}

double
mpc_map_get_tolerance()
{
    return 0;
}
