#include "rmap.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util/atomics.h"
#include "util/rmap.h"
#include "util/fast_hash.h"
#include "util/util.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256


// TODO: use fl_t to allow for using 128 bit floats (also requires d[2])
typedef union {
    double r;
    uint64_t d;
} bucket_t;

// float "equality" tolerance
static long double TOLERANCE = 1e-14l;
static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

struct rmap_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_float(bucket_t *b)
{
    printf("%.60f\n", b->r);
}

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("hex = %lx\n", b->d);
}

double
rmap_get_tolerance()
{
    return TOLERANCE;
}

static bool
real_close(double *in_table, const double* to_insert)
{
    return ((fabs(*in_table - *to_insert) < TOLERANCE));
}

int
rmap_find_or_put (const rmap_t *rmap, const double *v, ref_t *ret)
{
    bucket_t *val  = (bucket_t *) v;

    // Round the value to compute the hash with, but store the actual value v
    double round_v = roundl((long double) *v / TOLERANCE) * TOLERANCE;

    // fix 0 possibly having a sign
    if(round_v == 0.0) round_v = 0.0;
    
    uint32_t hash  = SuperFastHash(&round_v, sizeof(double), 0);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d != LOCK);
    assert (val->d != EMPTY);

    // Insert/lookup `v`
    for (unsigned int c = 0; c < rmap->threshold; c++) {
        size_t              ref = hash & rmap->mask;
        size_t              line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &rmap->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d == EMPTY) {
                if (cas(&bucket->d, EMPTY, LOCK)) {
                    *ret = ref;
                    atomic_write (&bucket->d, val->d);
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d) == LOCK) {}

            // 4. Bucket contains some complex value, check if close to `v`
            double *in_table = (double *)bucket;
            if (real_close(in_table, v)) {
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

double *
rmap_get (const rmap_t *rmap, const ref_t ref)
{
    return &rmap->table[ref].r;
}

uint64_t
rmap_count_entries (const rmap_t *rmap)
{
    uint64_t entries = 0;
    for (unsigned int c = 0; c < rmap->size; c++) {
        if (rmap->table[c].d != EMPTY)
            entries++;
    }
    return entries;
}

void
rmap_print_bitvalues(const rmap_t *rmap, const ref_t ref)
{
    bucket_t *b = (bucket_t *) rmap_get(rmap, ref);
    printf("%016lx", b->d);
}

rmap_t *
rmap_create (uint64_t size, double tolerance)
{
    TOLERANCE = tolerance;
    rmap_t  *rmap = calloc (1, sizeof(rmap_t));
    rmap->size = size;
    rmap->mask = rmap->size - 1;
    rmap->table = calloc (rmap->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < rmap->size; c++) {
        rmap->table[c].d = EMPTY;
    }
    rmap->threshold = rmap->size / 100;
    rmap->threshold = min(rmap->threshold, 1ULL << 16);
    rmap->seen_0 = 0;
    return rmap;
}

void
rmap_free (void *r)
{
    rmap_t * rmap = (rmap_t *) r;
    free (rmap->table);
    free (rmap);
}
