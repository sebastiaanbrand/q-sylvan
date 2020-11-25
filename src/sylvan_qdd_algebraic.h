#ifndef SYLVAN_QDD_ALGEBRAIC_H
#define SYLVAN_QDD_ALGEBRAIC_H

#endif

#include <stdint.h>

#include "sylvan_qdd_complex.h"

// 1 / sqrt(2)^k (a * omega^3 + b * omega^2 + c * omega + d)
typedef struct algebraic_s {
    int64_t a;
    int64_t b;
    int64_t c;
    int64_t d;
    int64_t k;
} algebraic_t;


algebraic_t algebraic_create(int64_t a, int64_t b, int64_t c, int64_t d, int64_t k);
algebraic_t algebraic_zero();
algebraic_t algebraic_one();
algebraic_t algebraic_add(algebraic_t a, algebraic_t b);
algebraic_t algebraic_mult(algebraic_t a, algebraic_t b);

void algebraic_init();
