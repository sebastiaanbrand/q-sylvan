#ifndef SYLVAN_QDD_COMPLEX_H
#define SYLVAN_QDD_COMPLEX_H

/**
 * Complex values with lookup table
 *
 * Based on [1]
 */

#include <math.h>
#include <stdint.h>

#include "util/cmap.h"
#include "util/rmap.h"
#include "util/tree_map.h"
#include "util/mpfr_tree_map.h"

typedef enum amp_storage_backend {
    COMP_HASHMAP, 
    REAL_HASHMAP, 
    REAL_TREE
} amp_storage_backend_t;

typedef uint64_t AMP;

AMP C_ZERO;
AMP C_ONE;
AMP C_MIN_ONE;

// GATE_ID's (gates are initialized in qdd_complex_init)
// currently 24 bits available for this number (see GATE_OPID)
static const uint32_t GATEID_I = 0;
static const uint32_t GATEID_X = 1;
static const uint32_t GATEID_Y = 2;
static const uint32_t GATEID_Z = 3;
static const uint32_t GATEID_H = 4;
static const uint32_t GATEID_S = 5;
static const uint32_t GATEID_T = 6;
static const uint32_t GATEID_Tdag = 7;
static const uint32_t GATEID_sqrtX = 8;
static const uint32_t GATEID_sqrtY = 9;
void init_gates();
// The next 255 gates are reserved for parameterized phase gates, which are
// initialized below. (these are GATEIDs 10 265)
// The reason why these are initialized beforhand instead of on-demand is that 
// we need a (for example) pi/16 gate always to have the same unique ID 
// throughout the entire run of the circuit.
// TODO: maybe find a better way?
void init_phase_gates(int n);
static inline uint32_t GATEID_Rk(int k) { return k + 10; };
// Another 255 parameterized phase gates, but this time with negative angles.
// (GATEIDs 266 through 521)
static inline uint32_t GATEID_Rk_dag(int k){ return k + 266; };

// Reserve 1000 gate IDs (comibned) for of the following: Rx, Ry, Rz. These 
// gate IDs are re-used when a user requires more than 1000 custom gates.
/**
 * Rotation around x-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Rx(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Rx(fl_t a);
/**
 * Rotation around y-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Ry(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Ry(fl_t a);
/**
 * Rotation around z-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Rz(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Rz(fl_t a);

static const uint64_t num_static_gates  = 522;
static const uint64_t num_dynamic_gates = 1000; // IDs in this rage are re-used 

// 2x2 gates, k := GATEID_U 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
AMP gates[522+1000][4]; // max 2^24 gates atm


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
AMP amp_div(AMP a, AMP b);

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

/* Comparing complex values */
bool comp_exact_equal(complex_t a, complex_t b);
bool comp_approx_equal(complex_t a, complex_t b);
bool comp_epsilon_close(complex_t a, complex_t b, long double epsilon);

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
uint64_t count_amplitude_table_enries();
uint64_t get_table_entries_estimate();
uint64_t get_table_size();
void free_amplitude_table();
void init_new_empty_table();
void delete_old_table();
AMP move_from_old_to_new(AMP);


#endif
