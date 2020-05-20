#include "cmap.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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
};

bool
cmap_find_or_put (const cmap_t *cmap, const complex_t *v, ref_t *ret)
{
    // TODO truncate

    //printf("r=%.60f,\ti=%.60f\n",v->r, v->i);
    
    // in some cases the following gives random values for some of the 256 bits
    // DONE: either switch to regular (64 bit) doubles, or identify these bits
    //bucket_t *b= (bucket_t *) v;
    //printf("value: %p %p\n", b->d[0], b->d[1]);
    /*
    complex_t trunc;
    trunc.r = v->r;
    trunc.i = v->i;

    bucket_t *buck = (bucket_t *) v;

    buck->d[0] = buck->d[0] & 0x0000000000000000;
    buck->d[1] = buck->d[1] & 0x0000000000000000;
    
    printf("value: \n%p\n%p\n", buck->d[0], buck->d[1]);
    
    printf("(%Lf, %Lf)\n",trunc.r, v->i);
    printf("15: %.15Lf\n",trunc.r, v->i);
    printf("20: %.20Lf\n",trunc.r, v->i);
    printf("25: %.25Lf\n",trunc.r, v->i);
    printf("30: %.30Lf\n",trunc.r, v->i);
    printf("35: %.35Lf\n",trunc.r, v->i);
    printf("40: %.40Lf\n",trunc.r, v->i);
    printf("45: %.45Lf\n",trunc.r, v->i);
    printf("50: %.50Lf\n",trunc.r, v->i);
    printf("60: %.60Lf\n",trunc.r, v->i);
    printf("65: %.65Lf\n",trunc.r, v->i);
    */

    uint32_t            hash = SuperFastHash(v, sizeof(complex_t), 0);
    uint32_t            prime = odd_primes[hash & PRIME_MASK];
    bucket_t *val = (bucket_t *) v;
    assert (val->d[0] != LOCK && val->d[0] != EMPTY);

    for (int c = 0; c < cmap->threshold; c++) {
        size_t              ref = hash & cmap->mask;
        size_t              line_end = (ref & CL_MASK) + CACHE_LINE_SIZE;
        for (size_t i = 0; i < CACHE_LINE_SIZE; i++) {
            bucket_t           *bucket = &cmap->table[ref];
            if (bucket->d[0] == EMPTY) {
                if (cas(&bucket->d[0], EMPTY, LOCK)) {
                    *ret = ref;
                    //atomic_write (&bucket->d[3], val->d[3]);
                    //atomic_write (&bucket->d[2], val->d[2]);
                    atomic_write (&bucket->d[1], val->d[1]);
                    atomic_write (&bucket->d[0], val->d[0]);
                    return 0;
                }
            }
            while (atomic_read(&bucket->d[0]) == LOCK) {}

            if (    bucket->d[0] == val->d[0] &&
                    bucket->d[1] == val->d[1] //&&
                    //bucket->d[2] == val->d[2] &&
                    //bucket->d[3] == val->d[3]
                    ) {
                *ret = ref;
                return 1;
            }
            ref += 1;
            ref = ref == line_end ? line_end - CACHE_LINE_SIZE : ref;
        }
        hash += prime << CACHE_LINE;
    }
    assert ("Hash table full" && false);
}

complex_t *
cmap_get (const cmap_t *cmap, const ref_t ref)
{
    return &cmap->table[ref].c;
}

cmap_t *
cmap_create (int size)
{
    cmap_t  *cmap = calloc (1, sizeof(cmap_t));
    cmap->size = 1ull << size;
    cmap->mask = cmap->size - 1;
    cmap->table = calloc (cmap->size, sizeof(bucket_t));
    for (int c = 0; c < cmap->size; c++) {
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
