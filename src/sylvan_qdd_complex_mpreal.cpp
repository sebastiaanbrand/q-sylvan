#include "sylvan_qdd_complex_mpreal.h"
#include "amp_storage/mpreal_tree_map.h"

#include "sylvan_int.h"
using namespace sylvan; // sylvan's internals are inside namespace in C++

// Table parameters
static const double default_tolerance = 1e-14;
static double tolerance;
size_t table_size;

// Actual table
mpreal_tree_map_t *mpreal_cmap;
mpreal_tree_map_t *mpreal_cmap_old;

static bool CACHE_AMP_OPS = true;

static mpfr::mpreal Pi;


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

bool mpreal_amp_exact_equal(AMP a, AMP b)
{
    return mpreal_comp_exact_equal(mpreal_comp_value(a), mpreal_comp_value(b));
}

bool mpreal_amp_approx_equal(AMP a, AMP b)
{
    return mpreal_comp_approx_equal(mpreal_comp_value(a), mpreal_comp_value(b));
}

bool mpreal_amp_epsilon_close(AMP a, AMP b, double epsilon)
{
    return mpreal_comp_epsilon_close(mpreal_comp_value(a), mpreal_comp_value(b), epsilon);
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

mpreal_complex
mpreal_comp_value_old(AMP a)
{
    // special cases
    if (a == C_ZERO)    return mpreal_comp_zero();
    if (a == C_ONE)     return mpreal_comp_one();
    if (a == C_MIN_ONE) return mpreal_comp_minus_one();

    // lookup
    mpreal_complex res;
    res = mpreal_tree_map_get(mpreal_cmap_old, a);
    return res;
}


AMP
mpreal_comp_lookup(mpreal_complex c)
{
    unsigned int res = 0;

    // TODO: local counters for gc instead of the counter in the mpreal_tree_map

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


mpreal_complex
mpreal_comp_add(mpreal_complex a, mpreal_complex b)
{
    mpreal_complex res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

mpreal_complex
mpreal_comp_sub(mpreal_complex a, mpreal_complex b)
{
    mpreal_complex res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

mpreal_complex
mpreal_comp_mul(mpreal_complex a, mpreal_complex b)
{
    mpreal_complex res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

mpreal_complex
mpreal_comp_div(mpreal_complex a, mpreal_complex b)
{
    mpreal_complex res;
    mpfr::mpreal denom;
    if (b.i == 0.0) {
        res.r = a.r / b.r;
        res.i = a.i / b.r;
    } else {
        denom = b.r * b.r + b.i * b.i;
        res.r = (a.r * b.r + a.i * b.i) / denom;
        res.i = (a.i * b.r - a.r * b.i) / denom;
    }
    return res;
}

double
mpreal_comp_to_prob(mpreal_complex a)
{
    mpfr::mpreal res = mpfr::sqrt ( (a.r*a.r) + (a.i*a.i) );
    res = res * res;
    return res.toDouble();
}

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

AMP
mpreal_amp_add(AMP a, AMP b)
{
    // special cases
    if (a == C_ZERO) return b;
    if (b == C_ZERO) return a;

    // check cache
    AMP res;
    if (CACHE_AMP_OPS) {
        if (cache_get3(CACHE_AMP_ADD, a, b, sylvan_false, &res)) {
            return res; // TODO: counters for these cache lookups/puts
        }
    }

    // compute and add result to ctable
    mpreal_complex ca, cb, cr;
    ca = mpreal_comp_value(a);
    cb = mpreal_comp_value(b);
    cr = mpreal_comp_add(ca, cb);
    res = mpreal_comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_ADD, a, b, sylvan_false, res);
    }
    return res;
}


AMP
mpreal_amp_sub(AMP a, AMP b)
{
    // special cases
    if (b == C_ZERO) return a;
    if (a == C_ZERO) return mpreal_amp_neg(b);

    // check cache
    AMP res;
    if (CACHE_AMP_OPS) {
        if (cache_get3(CACHE_AMP_SUB, a, b, sylvan_false, &res)) {
            return res; // TODO: counters for these cache lookups/puts
        }
    }

    // compute and add result to ctable
    mpreal_complex ca, cb, cr;
    ca = mpreal_comp_value(a);
    cb = mpreal_comp_value(b);
    cr = mpreal_comp_sub(ca, cb);
    res = mpreal_comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_SUB, a, b, sylvan_false, res);
    }
    return res;
}

AMP
mpreal_amp_mul(AMP a, AMP b)
{
    // special cases
    if (a == C_ONE) return b;
    if (b == C_ONE) return a;
    if (a == C_ZERO || b == C_ZERO) return C_ZERO;

    // check cache
    AMP res;
    if (CACHE_AMP_OPS) {
        if (cache_get3(CACHE_AMP_MUL, a, b, sylvan_false, &res)) {
            return res; // TODO: counters for these cache lookups/puts
        }
    }

    // compute and add result to ctable
    mpreal_complex ca, cb, cr;
    ca = mpreal_comp_value(a);
    cb = mpreal_comp_value(b);
    cr = mpreal_comp_mul(ca, cb);
    res = mpreal_comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_MUL, a, b, sylvan_false, res);
    }
    return res;
}

AMP
mpreal_amp_div(AMP a, AMP b)
{
    // special cases
    if (a == b)      return C_ONE;
    if (a == C_ZERO) return C_ZERO;
    if (b == C_ONE)  return a;

    // check cache
    AMP res;
    if (CACHE_AMP_OPS) {
        if (cache_get3(CACHE_AMP_DIV, a, b, sylvan_false, &res)) {
            return res; // TODO: counters for these cache lookups/puts
        }
    }

    // compute and hash result to ctable
    mpreal_complex ca, cb, cr;
    ca = mpreal_comp_value(a);
    cb = mpreal_comp_value(b);
    cr = mpreal_comp_div(ca, cb);
    res = mpreal_comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_DIV, a, b, sylvan_false, res);
    }
    return res;
}

double
mpreal_amp_to_prob(AMP a)
{
    return mpreal_comp_to_prob(mpreal_comp_value(a));
}

AMP
mpreal_prob_to_amp(double a)
{
    mpreal_complex c;
    c.r = mpfr::sqrt(a);
    c.i = 0;
    return mpreal_comp_lookup(c);
}


/* normalization of two amps */

AMP
mpreal_amp_normalize_low(AMP *low, AMP *high)
{
    // Normalize using low if low != 0
    AMP norm;
    if(*low != C_ZERO){
        mpreal_complex cl = mpreal_comp_value(*low);
        mpreal_complex ch = mpreal_comp_value(*high);
        ch    = mpreal_comp_div(ch, cl);
        *high = mpreal_comp_lookup(ch);
        norm  = *low;
        *low  = C_ONE;
    }
    else {
        norm  = *high;
        *high = C_ONE;
    }
    return norm;
}

AMP
mpreal_amp_normalize_largest(AMP *low, AMP *high)
{
    AMP norm;
    if (*low == *high) {
        norm  = *low;
        *low  = C_ONE;
        *high = C_ONE;
        return norm;
    }

    // Normalize using the absolute greatest value
    mpreal_complex cl = mpreal_comp_value(*low);
    mpreal_complex ch = mpreal_comp_value(*high);
    if ( (cl.r*cl.r + cl.i*cl.i)  >=  (ch.r*ch.r + ch.i*ch.i) ) {
        ch = mpreal_comp_div(ch, cl);
        *high = mpreal_comp_lookup(ch);
        norm = *low;
        *low  = C_ONE;
    }
    else {
        cl = mpreal_comp_div(cl, ch);
        *low = mpreal_comp_lookup(cl);
        norm  = *high;
        *high = C_ONE;
    }
    return norm;
}



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

    Pi = mpfr::const_pi();

    init_mpreal_gates();
}

void
init_mpreal_gates()
{
    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_ONE;

    k = GATEID_X;
    gates[k][0] = C_ZERO; gates[k][1] = C_ONE;
    gates[k][2] = C_ONE;  gates[k][3] = C_ZERO;

    k = GATEID_Y;
    gates[k][0] = C_ZERO; gates[k][1] = mpreal_comp_lookup(mpreal_comp_make(0.0, -1.0));
    gates[k][2] = mpreal_comp_lookup(mpreal_comp_make(0.0, 1.0));  gates[k][3] = C_ZERO;

    k = GATEID_Z;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = mpreal_comp_lookup(mpreal_comp_make(1.0/mpfr::sqrt(2),0));
    gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(-1.0/mpfr::sqrt(2),0));

    k = GATEID_S;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(0.0, 1.0));

    k = GATEID_T;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(1.0/mpfr::sqrt(2), 1.0/mpfr::sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(1.0/mpfr::sqrt(2), -1.0/mpfr::sqrt(2)));

    k = GATEID_sqrtX;
    gates[k][0] = mpreal_comp_lookup(mpreal_comp_make(0.5, 0.5)); gates[k][1] = mpreal_comp_lookup(mpreal_comp_make(0.5,-0.5));
    gates[k][2] = mpreal_comp_lookup(mpreal_comp_make(0.5,-0.5)); gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(0.5, 0.5));

    k = GATEID_sqrtY;
    gates[k][0] = mpreal_comp_lookup(mpreal_comp_make(0.5, 0.5)); gates[k][1] = mpreal_comp_lookup(mpreal_comp_make(-0.5,-0.5));
    gates[k][2] = mpreal_comp_lookup(mpreal_comp_make(0.5, 0.5)); gates[k][3] = mpreal_comp_lookup(mpreal_comp_make(0.5, 0.5));

    init_mpreal_phase_gates(255);

    //next_custom_id = 0;
}

void
init_mpreal_phase_gates(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    mpfr::mpreal angle;
    mpreal_complex cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / (1<<k);
        cartesian = mpreal_comp_make_angle(angle);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = mpreal_comp_lookup(cartesian);

        // backward rotation
        angle = -2*Pi / (double)(1<<k);
        cartesian = mpreal_comp_make_angle(angle);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = mpreal_comp_lookup(cartesian);
    }
}

void
free_mpreal_amplitude_table()
{
    mpreal_tree_map_free(mpreal_cmap);
}

void
init_new_empty_mpreal_table()
{
    // point old to current (full) ctable
    mpreal_cmap_old = mpreal_cmap;

    // re-init new (empty) ctable
    double tolerance = mpreal_tree_map_get_tolerance();
    init_mpreal_amplitude_table(table_size, tolerance);
}

void
delete_old_mpreal_table()
{
    // delete  old (full) table
    mpreal_tree_map_free(mpreal_cmap_old);
}

AMP
move_from_old_to_new_mpreal(AMP a)
{
    mpreal_complex c = mpreal_comp_value_old(a);
    return mpreal_comp_lookup(c);
}

uint64_t
get_mpreal_table_size()
{
    return table_size;
}

uint64_t
get_mpreal_table_num_entries()
{
    return mpreal_tree_map_get_num_entries(mpreal_cmap);
}

