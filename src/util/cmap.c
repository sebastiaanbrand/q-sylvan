#include "cmap.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util/atomics.h"
#include "util/cmap.h"
#include "util/fast_hash.h"
#include "util/util.h"

#undef CACHE_LINE
#undef CACHE_LINE_SIZE

#define CACHE_LINE 8
#define CACHE_LINE_SIZE 256


typedef union {
    complex_t       c;
    uint64_t        d[2];
} bucket_t;


static const uint64_t EMPTY = 14738995463583502973ull;
static const uint64_t LOCK  = 14738995463583502974ull;
static const uint64_t CL_MASK = -(1ULL << CACHE_LINE);

struct cmap_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};


void 
print_bucket_floats(bucket_t *b)
{
    printf("%.60f, %.60f\n", b->c.r, b->c.i);
}

void 
print_bucket_bits(bucket_t* b)
{
    printf("hex = %lx, %lx\n", b->d[0], b->d[1]);
}

bool
complex_close(complex_t *in_table, const complex_t* to_insert)
{
    return ((fabs(in_table->r - to_insert->r) < TOLERANCE) && 
            (fabs(in_table->i - to_insert->i) < TOLERANCE));
}

bool
cmap_find_or_put (const cmap_t *cmap, const complex_t *v, ref_t *ret)
{
    bucket_t *val  = (bucket_t *) v;

    // Round the value to compute the hash with, but store the actual value v
    complex_t round_v;
    round_v.r = round(v->r * INV_TOLERANCE)/INV_TOLERANCE;
    round_v.i = round(v->i * INV_TOLERANCE)/INV_TOLERANCE;

    // fix 0 possibly having a sign
    if(round_v.r == 0.0) round_v.r = 0.0;
    if(round_v.i == 0.0) round_v.i = 0.0;
    
    uint32_t hash  = SuperFastHash(&round_v, sizeof(complex_t), 0);
    uint32_t prime = odd_primes[hash & PRIME_MASK];

    assert (val->d[0] != LOCK);
    assert (val->d[0] != EMPTY);

    // Insert/lookup `v`
    for (unsigned int c = 0; c < cmap->threshold; c++) {
        size_t              ref = hash & cmap->mask;
        size_t              line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            
            // 1. Get bucket
            bucket_t *bucket = &cmap->table[ref];

            // 2. If bucket empty, insert new value here
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    atomic_write (&bucket->d[1], val->d[1]);
                    atomic_write (&bucket->d[0], val->d[0]);
                    return 0;
                }
            }

            // 3. Bucket not empty, wait for lock
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            // 4. Bucket contains some complex value, check if close to `v`
            complex_t *in_table = (complex_t *)bucket;
            if (complex_close(in_table, v)) {
                *ret = ref;
                return 1;
            }

            // If unsuccessful, try next
            ref += 1;
            ref = ref == line_end ? line_end - CACHE_LINE_SIZE : ref;
        }
        hash += prime << CACHE_LINE;
    }
    assert ("Hash table full" && false);
    return -1; // (should be) unreachable but deals with compiler warning
}

complex_t *
cmap_get (const cmap_t *cmap, const ref_t ref)
{
    return &cmap->table[ref].c;
}

void
print_bitvalues(const cmap_t *cmap, const ref_t ref)
{
    bucket_t *b = (bucket_t *) cmap_get(cmap, ref);
    printf("%016lx %016lx", b->d[0], b->d[1]);
}

cmap_t *
cmap_create (int size)
{
    cmap_t  *cmap = calloc (1, sizeof(cmap_t));
    cmap->size = 1ull << size;
    cmap->mask = cmap->size - 1;
    cmap->table = calloc (cmap->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < cmap->size; c++) {
        cmap->table[c].d[0] = EMPTY;
    }
    cmap->threshold = cmap->size / 100;
    cmap->threshold = min(cmap->threshold, 1ULL << 16);
    cmap->seen_0 = 0;
    return cmap;
}

void
cmap_free (cmap_t *cmap)
{
    free (cmap->table);
    free (cmap);
}
