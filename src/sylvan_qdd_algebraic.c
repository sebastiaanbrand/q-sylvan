/**
 * Algebraic encoding of amplitudes, as in 
 * "Overcoming the Tade-off between Accuracy and Compactness in Decision 
 * Diagrams for Quantum Computation", P. Niemann, R. Drechsler (2020)
 * 
 */

#include <stdio.h>

#include "sylvan_qdd_algebraic.h"


static long double Pi;


complex_t omega;
complex_t omega_2;
complex_t omega_3;

algebraic_t
algebraic_create(int64_t a, int64_t b, int64_t c, int64_t d, int64_t k)
{
    algebraic_t res;
    res.a = a;
    res.b = b;
    res.c = c;
    res.d = d;
    res.k = k;
    return res;
}



algebraic_t
algebraic_zero()
{
    return algebraic_create(0, 0, 0, 0, 0);
}

algebraic_t
algebraic_one() 
{
    return algebraic_create(0, 0, 0, 1, 0);
}

complex_t
algebraic_to_comp(algebraic_t a)
{
    complex_t res = comp_zero();
    res = comp_add(res, comp_mul(omega_3, comp_make((double)a.a, 0.0)));
    res = comp_add(res, comp_mul(omega_2, comp_make((double)a.b, 0.0)));
    res = comp_add(res, comp_mul(omega,   comp_make((double)a.c, 0.0)));
    res = comp_add(res, comp_make((double)a.d, 0.0));
    res = comp_mul(res, comp_make(  1.0 / (pow(sqrt(2.0),(double)a.k)) , 0.0));
    return res;
}

algebraic_t
algebraic_add(algebraic_t a, algebraic_t b)
{
    // TODO
    return algebraic_zero();
}

algebraic_t
algebraic_mult(algebraic_t a, algebraic_t b)
{
    // TODO
    return algebraic_zero();
}

void
algebraic_init()
{
    Pi = 2.0 * acos(0.0);

    omega   = comp_make(1.0/sqrt(2.0), 1.0/sqrt(2.0));
    omega_2 = comp_make(0.0, 1.0);
    omega_3 = comp_make(-1.0/sqrt(2.0), 1.0/sqrt(2.0));
}
