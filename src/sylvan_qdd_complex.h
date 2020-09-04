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


//typedef struct complex_s {
//   long double r;   // real
//   long double i;   // imaginary
//} complex_t;


typedef uint64_t AMP;

AMP C_ZERO;
AMP C_ONE;

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

// 2x2 gates, k := GATEID_U 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
AMP gates[522][4]; // max 2^16 gates atm


//uint32_t Ctentries;

// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

// Some of these methods are exposed for unit testing, but don't need to be
// used in the qdd implementation.
complex_t Cmake (double r, double i);
double Qmake (int a, int b, int c);
AMP Clookup (complex_t c);

complex_t Cvalue(AMP);
void Cprint(complex_t);
void Cprint_bitvalues(AMP);

bool comp_exact_equal(complex_t a, complex_t b);
bool comp_approx_equal(complex_t a, complex_t b);
bool comp_epsilon_close(complex_t a, complex_t b, double epsilon);


/* Arithmetic operations on AMPs */
AMP amp_neg(AMP a);
AMP amp_add(AMP a, AMP b);
AMP amp_sub(AMP a, AMP b);
AMP amp_mul(AMP a, AMP b);
AMP amp_div(AMP a, AMP b);
//cint CAbs(cint); /// by PN: returns the absolut value of a complex number
//cint CUnit(cint a); ///by PN: returns whether a complex number has norm 1

/* Arithmetic operations on complex structs */
complex_t comp_negative(complex_t a);

//bool Ccomp(complex_t x, complex_t y);

void init_amplitude_table(size_t size);
uint64_t count_amplitude_table_enries();
void free_amplitude_table();
void init_new_empty_table();
void delete_old_table();
AMP move_from_old_to_new(AMP);


#endif
