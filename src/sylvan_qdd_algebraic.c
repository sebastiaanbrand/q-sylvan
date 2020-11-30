/**
 * Algebraic encoding of amplitudes, as in 
 * "Overcoming the Tade-off between Accuracy and Compactness in Decision 
 * Diagrams for Quantum Computation", P. Niemann, R. Drechsler (2020)
 * 
 */

#include <stdio.h>
#include <assert.h>

#include "sylvan_qdd_algebraic.h"


complex_t omega;
complex_t omega_2;
complex_t omega_3;


algebraic_t
algebraic_pack(int64_t k, int64_t a, int64_t b, int64_t c, int64_t d)
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
algebraic_create(int64_t k, int64_t a, int64_t b, int64_t c, int64_t d)
{
    algebraic_t res = algebraic_pack(k, a, b, c, d);
    return algebraic_minimal(res);
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


algebraic_t
algebraic_move_sqrt2_inside(algebraic_t x)
{
    algebraic_t res;

    // take factor sqrt(2) out of first term
    res.k = x.k + 1;

    // multiply (a,b,c,d) with sqrt(2) = (a=-1, b=0, c=1, d=0)
    res.a = x.b - x.d;
    res.b = x.c + x.a;
    res.c = x.b + x.d;
    res.d = x.c - x.a;

    return res;
}


algebraic_t
algebraic_add(algebraic_t x, algebraic_t y)
{
    // TODO: LaTeX mathematical writeup of this function.
    algebraic_t res, temp_x, temp_y;

    if ( (x.k - y.k) % 2 != 0) {
        // either x.k or y.k is odd. Out of x and y, multiply the one 
        // with the odd exponent with sqrt(2) by changing (a,b,c,d),
        // rather than k, such that we can increase k by 1 and make it even.
        if (x.k % 2 != 0) {
            assert(y.k % 2 == 0);
            temp_x = algebraic_move_sqrt2_inside(x);
            temp_y = y;
        }
        else {
            assert(y.k % 2 != 0);
            temp_x = x;
            temp_y = algebraic_move_sqrt2_inside(y);
        }
    }
    else {
        temp_x = x;
        temp_y = y;
    }
    x = temp_x;
    y = temp_y;

    assert((x.k - y.k) % 2 == 0);

    // If (x.k) and (y.k) differ by a multiple of 2, the first factor of 
    // x and y differ from each other by a factor which is an integer 
    // power of 2.
    // The following orders x and y such that this 2^z is an integer.
    if (x.k < y.k) {
        temp_x = y;
        temp_y = x;
    }
    x = temp_x;
    y = temp_y;

    assert((x.k - y.k) >= 0); // x.k - y.k positive and even so norm = int
    uint64_t norm = 1<<((x.k - y.k)/2);
    // the following takes the factor 'norm' inside the backets for y,
    // such that x and y will both have the same factor up front
    y.k = x.k;
    y.a = y.a * norm;
    y.b = y.b * norm;
    y.c = y.c * norm;
    y.d = y.d * norm;

    // now the factor in front is the same, 
    // so we can just add them linearly
    assert(x.k == y.k);
    res.k = (x.k);
    res.a = (x.a + y.a);
    res.b = (x.b + y.b);
    res.c = (x.c + y.c);
    res.d = (x.d + y.d);
    
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
