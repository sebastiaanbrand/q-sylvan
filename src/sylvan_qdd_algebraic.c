/**
 * Algebraic encoding of amplitudes, as in 
 * "Overcoming the Tade-off between Accuracy and Compactness in Decision 
 * Diagrams for Quantum Computation", P. Niemann, R. Drechsler (2020)
 * 
 */

#include <stdio.h>

#include "sylvan_qdd_algebraic.h"


complex_t omega;
complex_t omega_2;
complex_t omega_3;

algebraic_t
algebraic_create(int64_t k, int64_t a, int64_t b, int64_t c, int64_t d)
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
    return algebraic_create(0, 0, 0, 0, 1);
}


algebraic_t
algebraic_sqrt2(int64_t k)
{
    return algebraic_create(-k, 0, 0, 0, 1);
}


algebraic_t
algebraic_minimal(algebraic_t x)
{
    // TODO
    return x;
}



algebraic_t
algebraic_add(algebraic_t x, algebraic_t y)
{
    // TODO
    return algebraic_zero();
}


algebraic_t
algebraic_mult(algebraic_t x, algebraic_t y)
{
    algebraic_t res;
    res.k =  (x.k + y.k);
    res.a =  (x.a * y.d) + (x.b * y.c) + (x.c * y.b) + (x.d * y.a);
    res.b = -(x.a * y.a) + (x.b * y.d) + (x.c * y.c) + (x.d * y.b);
    res.c = -(x.a * y.b) - (x.b * y.a) + (x.c * y.d) + (x.d * y.c);
    res.d = -(x.a * y.c) - (x.b * y.b) - (x.c * y.a) + (x.d * y.d);
    return algebraic_minimal(res);
}


complex_t
algebraic_to_comp(algebraic_t x)
{
    complex_t res = comp_zero();
    res = comp_add(res, comp_mul(omega_3, comp_make((double)x.a, 0.0)));
    res = comp_add(res, comp_mul(omega_2, comp_make((double)x.b, 0.0)));
    res = comp_add(res, comp_mul(omega,   comp_make((double)x.c, 0.0)));
    res = comp_add(res, comp_make((double)x.d, 0.0));
    res = comp_mul(res, comp_make(  1.0 / (pow(sqrt(2.0),(double)x.k)) , 0.0));
    return res;
}


void
algebraic_init()
{
    omega   = comp_make(1.0/sqrt(2.0), 1.0/sqrt(2.0));
    omega_2 = comp_make(0.0, 1.0);
    omega_3 = comp_make(-1.0/sqrt(2.0), 1.0/sqrt(2.0));
}
