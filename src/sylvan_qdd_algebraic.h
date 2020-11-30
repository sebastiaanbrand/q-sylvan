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


algebraic_t algebraic_create(int64_t k, int64_t a, int64_t b, int64_t c, int64_t d);

/* Algebraic representation of 0 */
algebraic_t algebraic_zero();

/* Algebraic representation of 1 */
algebraic_t algebraic_one();

/* Algebraic representation of sqrt(2)^k */
algebraic_t algebraic_sqrt2(int64_t k);

algebraic_t algebraic_minimal(algebraic_t x);

/* Returns the minimized algebraic representation of x + y */
algebraic_t algebraic_add(algebraic_t x, algebraic_t y);

/* Returns the minimized algebraic representation of x * y */
algebraic_t algebraic_mult(algebraic_t x, algebraic_t y);

/* Returns the complex number encoded in x in cartesian coordinates form */
complex_t algebraic_to_comp(algebraic_t x);

void algebraic_init();
