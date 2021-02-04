#include "sylvan_qdd_complex_mpreal.h"
#include "util/mpreal_tree_map.h"

// Table parameters
static const double default_tolerance = 1e-14;
static double tolerance;
size_t table_size;

// Actual table
mpreal_tree_map_t *mpreal_cmap;
mpreal_tree_map_t *mpreal_cmap_old;


/* Comparing complex values */

bool mpreal_comp_exact_equal(mpreal_complex a, mpreal_complex b)
{
    return (a.r == b.r && a.i == b.i);
}

bool mpreal_comp_approx_equal(mpreal_complex a, mpreal_complex b)
{
    return mpreal_comp_epsilon_close(a, b, mpreal_tree_map_get_tolerance());
}

bool mpreal_comp_epsilon_close(mpreal_complex a, mpreal_complex b, double epsilon)
{
    return ( (mpfr::abs(a.r - b.r) < epsilon) && (mpfr::abs(a.i - b.i) < epsilon) );
}


/* Shorthand functions for making complex numbers */

mpreal_complex
mpreal_comp_make(mpfr::mpreal r, mpfr::mpreal i)
{
    mpreal_complex res;
    res.r = r;
    res.i = i;
    return res;
}


mpreal_complex
mpreal_comp_make_angle(mpfr::mpreal theta)
{
    mpreal_complex res;
    res.r = mpfr::cos(theta);
    res.i = mpfr::sin(theta);
    return res;
}

mpreal_complex
mpreal_comp_zero()
{
    return mpreal_comp_make(0.0, 0.0);
}

mpreal_complex
mpreal_comp_one()
{
    return mpreal_comp_make(1.0, 0.0);
}

mpreal_complex
mpreal_comp_minus_one()
{
    return mpreal_comp_make(-1.0, 0.0);
}


mpfr::mpreal mpreal_sqrt2(mpfr::mpreal a)
{
    mpfr::mpreal aa = a;
    return aa * mpfr::sqrt(2);
}

/* value (get) cand lookup (find or put) */

mpreal_complex
mpreal_comp_value(AMP a)
{
    // special cases
    if (a == C_ZERO)    return mpreal_comp_zero();
    if (a == C_ONE)     return mpreal_comp_one();
    if (a == C_MIN_ONE) return mpreal_comp_minus_one();

    // lookup
    mpreal_complex res;
    res = mpreal_tree_map_get(mpreal_cmap, a);
    return res;
}

AMP
mpreal_comp_lookup(mpreal_complex c)
{
    unsigned int res = 0;

    // TODO: gc buffer counters

    mpreal_tree_map_find_or_put(mpreal_cmap, c, &res);

    return res;
}


/* Arithmetic operations on mpreal complex structs */

mpreal_complex
mpreal_comp_abs(mpreal_complex a)
{
    mpreal_complex res;
    res.r = mpfr::sqrt( (a.r*a.r)  + (a.i*a.i) );
    res.i = 0;
    return res;
}


mpreal_complex
mpreal_comp_neg(mpreal_complex a)
{
    mpreal_complex res;
    res.r = -a.r;
    res.i = -a.i;
    return res;
}

/*
mpreal_complex
mpreal_comp_add(mpreal_complex a, mpreal_complex b)
{
    // TODO
    return 0;
}

mpreal_complex
mpreal_comp_sub(mpreal_complex a, mpreal_complex b)
{
    // TODO
    return 0;
}

mpreal_complex
mpreal_comp_mul(mpreal_complex a, mpreal_complex b)
{
    // TODO
    return 0;
}

mpreal_complex
mpreal_comp_div(mpreal_complex a, mpreal_complex b)
{
    // TODO
    return 0;
}
*/


/* Arithmetic operations on AMPs */
AMP
mpreal_amp_abs(AMP a)
{
     // special cases
    if (a == C_ZERO || a == C_ONE) return a;
    if (a == C_MIN_ONE) return C_ONE;

    mpreal_complex ca, cr;
    AMP res;

    ca = mpreal_comp_value(a);
    cr = mpreal_comp_abs(ca);

    res = mpreal_comp_lookup(cr);
    return res;
}


AMP
mpreal_amp_neg(AMP a)
{
    // special cases
    if (a == C_ZERO) return C_ZERO;
    if (a == C_ONE) return C_MIN_ONE;
    if (a == C_MIN_ONE) return C_ONE;

    mpreal_complex ca, cr;
    AMP res;

    ca = mpreal_comp_value(a);
    cr = mpreal_comp_neg(ca);

    res = mpreal_comp_lookup(cr);
    return res;
}

/*

AMP
mpreal_amp_add(AMP a, AMP b)
{
    // TODO
    return 0;
}

AMP
mpreal_amp_sub(AMP a, AMP b)
{
    // TODO
    return 0;
}

AMP
mpreal_amp_mul(AMP a, AMP b)
{
    // TODO
    return 0;
}

AMP
mpreal_amp_div(AMP a, AMP b)
{
    // TODO
    return 0;
}
*/



/* Managing the complex value table */

void
init_mpreal_amplitude_table(size_t size, long double tol)
{
    tolerance = (tol < 0) ? default_tolerance : tol;
    table_size = size;

    //table_entries_est = 0;
    //table_entries_local = 0;

    mpreal_cmap = mpreal_tree_map_create(table_size, tolerance);
    
    // NOTE: the sum of the local counters sometimes exceeds the actual total
    // number of entries (when just counting the global value with atomic adds
    // after every insert). This might be because 'ctable_entries_local = 0' 
    // doesn't set it to 0 for all threads. Since it's only an estimate it is
    // not a huge issue though.
    // TODO: figure out how to get lace to handle this counting better


    C_ONE     = mpreal_comp_lookup(mpreal_comp_one());
    C_ZERO    = mpreal_comp_lookup(mpreal_comp_zero());
    C_MIN_ONE = mpreal_comp_lookup(mpreal_comp_minus_one());

    //Pi = 2.0 * flt_acos(0.0);

    //init_gates();
}

void
free_mpreal_amplitude_table()
{
    mpreal_tree_map_free(mpreal_cmap);
}