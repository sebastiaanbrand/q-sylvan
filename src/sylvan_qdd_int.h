#ifndef SYLVAN_QDD_INT_H
#define SYLVAN_QDD_INT_H

/**
 * Complex values with lookup table
 *
 * Based on [1]
 */


#include <math.h>
#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct complex_s {
   long double r;   // real
   long double i;   // imaginary
   long double m;   // magnitude
   long double a;   // angle
} complex_t;


typedef uint32_t cint;

static const cint C_ZERO = 0; // not hashed
static const cint C_ONE  = 1; // not hashed

// 2x2 gates, 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
cint gates[10][4];

uint32_t Ctentries;

// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

// Some of these methods are exposed for unit testing, but don't need to be
// used in the qdd implementation.
complex_t Cmake (long double r, long double i);
cint Clookup (complex_t c);

complex_t Cvalue(cint);
void Cprint(complex_t);

cint Cnegative(cint);
cint Cadd(cint,cint);
cint Csub(cint,cint);
cint Cmul(cint,cint);
cint CintMul(cint,cint); // multiply by an integer
cint Cdiv(cint,cint);
cint CAbs(cint); /// by PN: returns the absolut value of a complex number
cint CUnit(cint a); ///by PN: returns whether a complex number has norm 1

bool Ccomp(complex_t x, complex_t y);



// TODO: call somewhere
void qdd_complex_init();

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
