#ifndef SYLVAN_QDD_COMPLEX_MPREAL_H
#define SYLVAN_QDD_COMPLEX_MPREAL_H

#include <stdint.h>

#include "util/mpfr_tree_map.h"

typedef uint64_t AMP;

#ifdef __cplusplus
extern "C" {
#endif

int mpreal_test();

/* Arithmetic operations on AMPs */
AMP mpreal_amp_abs(AMP a);
AMP mpreal_amp_neg(AMP a);
AMP mpreal_amp_add(AMP a, AMP b);
AMP mpreal_amp_sub(AMP a, AMP b);
AMP mpreal_amp_mul(AMP a, AMP b);
AMP mpreal_amp_div(AMP a, AMP b);

/* Managing the complex value table */
void init_mpreal_amplitude_table(size_t size, long double tol);
void free_mpreal_amplitude_table();

#ifdef __cplusplus
}
#endif

#endif // SYLVAN_QDD_COMPLEX_MPREAL_H
