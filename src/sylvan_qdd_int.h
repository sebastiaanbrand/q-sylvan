#ifndef SYLVAN_QDD_INT_H
#define SYLVAN_QDD_INT_H

/**
 * Complex values with lookup table
 *
 * Based on [1]
 */


#include <math.h>
#include <stdint.h>


typedef uint32_t cint;

static const cint CZRO = 0; // not hashed
static const cint CONE = 1; // not hashed

// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

cint Cnegative(cint);
cint Cadd(cint,cint);
cint Csub(cint,cint);
cint Cmul(cint,cint);
cint CintMul(cint,cint); // multiply by an integer
cint Cdiv(cint,cint);
cint CAbs(cint); /// by PN: returns the absolut value of a complex number
cint CUnit(cint a); ///by PN: returns whether a complex number has norm 1

// TODO: call somewhere
void qdd_complex_init();

#endif
