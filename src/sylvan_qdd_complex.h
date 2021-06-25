#ifndef SYLVAN_QDD_COMPLEX_H
#define SYLVAN_QDD_COMPLEX_H

/**
 * Complex values with lookup table
 *
 * Based on [1]
 */

#include <math.h>
#include <stdint.h>

#include "amp_storage/amp_storage_interface.h"


typedef uint64_t AMP;

extern const bool CACHE_AMP_OPS;
extern const bool CACHE_INV_OPS;

AMP C_ZERO;
AMP C_ONE;
AMP C_MIN_ONE;


/* Shorthand functions for making complex numbers */
complex_t comp_make(fl_t r, fl_t i);
complex_t comp_make_angle(fl_t theta);
complex_t comp_zero();
complex_t comp_one();
complex_t comp_minus_one();
fl_t comp_qmake(int a, int b, int c);

/* Arithmetic operations on AMPs */
AMP amp_abs(AMP a);
AMP amp_neg(AMP a);
AMP amp_add(AMP a, AMP b);
AMP amp_sub(AMP a, AMP b);
AMP amp_mul(AMP a, AMP b);
AMP amp_mul_down(AMP a, AMP b); // same as mul but uses different stats counter for propagating edge weights down
AMP amp_div(AMP a, AMP b);
double amp_to_prob(AMP a);
AMP prob_to_amp(double a); // AMP of (sqrt(a), 0)

/* Arithmetic operations on complex structs */
complex_t comp_abs(complex_t a);
complex_t comp_neg(complex_t a);
complex_t comp_cnj(complex_t a);
complex_t comp_add(complex_t a, complex_t b);
complex_t comp_sub(complex_t a, complex_t b);
complex_t comp_mul(complex_t a, complex_t b);
complex_t comp_div(complex_t a, complex_t b);
double comp_to_prob(complex_t a);

/* normalization of two amps */
AMP amp_normalize_low(AMP *low, AMP *high);
AMP amp_normalize_largest(AMP *low, AMP *high);
#define amp_normalize_low_ptr &amp_normalize_low
#define amp_normalize_largest_ptr &amp_normalize_largest

/* Comparing AMPs */
bool amp_exact_equal(AMP a, AMP b);
bool amp_approx_equal(AMP a, AMP b);
bool amp_epsilon_close(AMP a, AMP b, double epsilon);

/* Comparing complex values */
bool comp_exact_equal(complex_t a, complex_t b);
bool comp_approx_equal(complex_t a, complex_t b);
bool comp_epsilon_close(complex_t a, complex_t b, double epsilon);

/* Inserting / retrieving complex values from complex table */
AMP comp_lookup(complex_t c);
AMP comp_try_lookup(complex_t c, bool *success);
complex_t comp_value(AMP a);

/* Printing */
void comp_print(complex_t c);
void comp_print_sci(complex_t c);
void comp_print_digits(complex_t c, uint32_t digits);
void comp_print_digits_sci(complex_t c, uint32_t digits);
void comp_print_bits(AMP a);

/* Managing the complex value table */
void init_amplitude_table(size_t size, long double tolerance, amp_storage_backend_t backend);
double amp_store_get_tolerance();
uint64_t count_amplitude_table_enries();
uint64_t get_table_entries_estimate();
uint64_t get_table_size();
void free_amplitude_table();
void init_new_empty_table();
void delete_old_table();
AMP move_from_old_to_new(AMP);


#endif
