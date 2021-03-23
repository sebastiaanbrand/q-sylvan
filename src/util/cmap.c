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

// how many "blocks" of 64 bits for a single table entry
#define entry_size (2*sizeof(fl_t)/8)

typedef union {
    complex_t       c;
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
typedef struct cmap_s cmap_t;
struct cmap_s {
    size_t              size;
    size_t              mask;
    size_t              threshold;
    int                 seen_0;
    bucket_t  __attribute__(( __aligned__(32)))       *table;
    // Q: should this 32 change to 16 now that we use doubles instead of
    // long doubles for the real and imaginary components?
};

static void __attribute__((unused))
print_bucket_floats(bucket_t *b)
{
    printf("%.60Lf, %.60Lf\n", (long double) b->c.r, (long double) b->c.i);
}

static void __attribute__((unused))
print_bucket_bits(bucket_t* b)
{
    printf("%016lx", b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016lx", b->d[k]);
    }
    printf("\n");
}

double
cmap_get_tolerance()
{
    return TOLERANCE;
}

static bool
complex_close(complex_t *in_table, const complex_t* to_insert)
{
    return ((flt_abs(in_table->r - to_insert->r) < TOLERANCE) && 
            (flt_abs(in_table->i - to_insert->i) < TOLERANCE));
}

int
cmap_find_or_put(const void *dbs, const complex_t *v, ref_t *ret)
{
    cmap_t *cmap = (cmap_t *) dbs;
    bucket_t *val  = (bucket_t *) v;

    // Round the value to compute the hash with, but store the actual value v
    bucket_t round_v;
    round_v.c.r = flt_round(v->r / TOLERANCE) * TOLERANCE;
    round_v.c.i = flt_round(v->i / TOLERANCE) * TOLERANCE;

    // fix 0 possibly having a sign
    if(round_v.c.r == 0.0) round_v.c.r = 0.0;
    if(round_v.c.i == 0.0) round_v.c.i = 0.0;

    //printf("(%.3f,%.3f) ",(float)round_v.c.r,(float)round_v.c.i);
    //print_bucket_bits(&round_v); 
    
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
    // amplitude table full, unable to add
    return -1;
}

complex_t *
cmap_get(const void *dbs, const ref_t ref)
{
    cmap_t *cmap = (cmap_t *) dbs;
    return &cmap->table[ref].c;
}

uint64_t
cmap_count_entries(const void *dbs)
{
    cmap_t *cmap = (cmap_t *) dbs;
    uint64_t entries = 0;
    for (unsigned int c = 0; c < cmap->size; c++) {
        if (cmap->table[c].d[0] != EMPTY)
            entries++;
    }
    return entries;
}

void
print_bitvalues(const void *dbs, const ref_t ref)
{
    cmap_t *cmap = (cmap_t *) dbs;
    bucket_t *b = (bucket_t *) cmap_get(cmap, ref);
    printf("%016lx", b->d[0]);
    for (unsigned int k = 1; k < entry_size; k++) {
        printf(" %016lx", b->d[k]);
    }
}

void *
cmap_create(uint64_t size, double tolerance)
{
    TOLERANCE = tolerance;
    cmap_t  *cmap = calloc (1, sizeof(cmap_t));
    cmap->size = size;
    cmap->mask = cmap->size - 1;
    cmap->table = calloc (cmap->size, sizeof(bucket_t));
    for (unsigned int c = 0; c < cmap->size; c++) {
        cmap->table[c].d[0] = EMPTY;
    }
    cmap->threshold = cmap->size / 100;
    cmap->threshold = min(cmap->threshold, 1ULL << 16);
    cmap->seen_0 = 0;
    return (void *) cmap;
}

void
cmap_free(void *dbs)
{
    cmap_t * cmap = (cmap_t *) dbs;
    free (cmap->table);
    free (cmap);
}
