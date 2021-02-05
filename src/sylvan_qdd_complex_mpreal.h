#ifndef SYLVAN_QDD_COMPLEX_MPREAL_H
#define SYLVAN_QDD_COMPLEX_MPREAL_H

#include <stdint.h>

#include "sylvan_qdd_complex.h" // using for C_ONE/C_ZERO (ugly but we'll clean up later)


typedef uint64_t AMP;


// C++ exclusive functions
#ifdef __cplusplus 
#include "util/mpreal_tree_map.h"

/* Comparing complex values */
bool mpreal_comp_exact_equal(mpreal_complex a, mpreal_complex b);
bool mpreal_comp_approx_equal(mpreal_complex a, mpreal_complex b);
bool mpreal_comp_epsilon_close(mpreal_complex a, mpreal_complex b, double epsilon);

/* value (get) cand lookup (find or put) */
mpreal_complex mpreal_comp_value(AMP a); // get
AMP mpreal_comp_lookup(mpreal_complex c); // find or put

/* Shorthand functions for making complex numbers */
//mpreal_complex mpreal_comp_make(double r, double i);
mpreal_complex mpreal_comp_make(mpfr::mpreal r, mpfr::mpreal i);
mpreal_complex mpreal_comp_make_angle(double theta);
mpreal_complex mpreal_comp_zero();
mpreal_complex mpreal_comp_one();
mpreal_complex mpreal_comp_minus_one();
mpfr::mpreal mpreal_sqrt2(mpfr::mpreal a); // a * sqrt2

/* Arithmetic operations on mpreal complex structs */
mpreal_complex mpreal_comp_abs(mpreal_complex a);
mpreal_complex mpreal_comp_neg(mpreal_complex a);
mpreal_complex mpreal_comp_add(mpreal_complex a, mpreal_complex b);
mpreal_complex mpreal_comp_sub(mpreal_complex a, mpreal_complex b);
mpreal_complex mpreal_comp_mul(mpreal_complex a, mpreal_complex b);
mpreal_complex mpreal_comp_div(mpreal_complex a, mpreal_complex b);


#endif



// Fuctions which can be called from a C file.
#ifdef __cplusplus
extern "C" {
#endif

/* Arithmetic operations on AMPs */
AMP mpreal_amp_abs(AMP a);
AMP mpreal_amp_neg(AMP a);
AMP mpreal_amp_add(AMP a, AMP b);
AMP mpreal_amp_sub(AMP a, AMP b);
AMP mpreal_amp_mul(AMP a, AMP b);
AMP mpreal_amp_div(AMP a, AMP b);

/* normalization of two amps */
AMP mpreal_amp_normalize_low(AMP *low, AMP *high);
AMP mpreal_amp_normalize_largest(AMP *low, AMP *high);

/* Managing the complex value table */
void init_mpreal_amplitude_table(uint64_t size, long double tol);
void free_mpreal_amplitude_table();
void init_mpreal_gates();
void init_mpreal_phase_gates(int n);

#ifdef __cplusplus
}
#endif


#endif // SYLVAN_QDD_COMPLEX_MPREAL_H
