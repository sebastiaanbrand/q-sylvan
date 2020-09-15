/***************************************************************

Adapted from implementation by:

January 28, 2008
Michael Miller
University of Victoria
Victoria, BC 
CANADA V8W 3P6
mmiller@cs.uvic.ca

****************************************************************/

#include <stdio.h>

#include "sylvan_int.h"
#include "sylvan_qdd_complex.h"
#include "util/cmap.h"


/**********************************************
Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation
**********************************************/
#define qdd_cos(fac,div) cosl((long double)(fac)*Pi/(long double)(div))
#define qdd_sin(fac,div) sinl((long double)(fac)*Pi/(long double)(div))


static long double Pi;    // set value of global Pi


size_t SIZE;
cmap_t *ctable;
cmap_t *ctable_old;
static bool CACHE_AMP_OPS = true;



/* Shorthand functions for making complex numbers */

complex_t
comp_make(double r, double i)
{
    complex_t res;
    res.r = r;
    res.i = i;
    return res;
}

complex_t
comp_make_angle(double theta)
{
    complex_t c;
    c.r = cos(theta);
    c.i = sin(theta);
    return c;
}

complex_t
comp_zero()
{
    return comp_make(0.0, 0.0);
}

complex_t
comp_one()
{
    return comp_make(1.0, 0.0);
}

complex_t
comp_minus_one()
{
    return comp_make(-1.0, 0.0);
}



double
comp_qmake(int a, int b, int c)
{
    return (((double) a + ((double) b) * sqrt (2.0)) / (double) (c));
}



/* Arithmetic operations on AMPs */

AMP
amp_abs(AMP a)
{
    // special cases
    if (a == C_ZERO || a == C_ONE) return a;
    if (a == C_MIN_ONE) return C_ONE;

    complex_t ca, cr;
    AMP res;

    ca = comp_value(a);
    cr = comp_abs(ca);

    res = comp_lookup(cr);
    return res;
}

AMP
amp_neg(AMP a)
{
    // special cases
    if (a == C_ZERO) return C_ZERO;
    if (a == C_ONE) return C_MIN_ONE;
    if (a == C_MIN_ONE) return C_ONE;

    complex_t ca, cr;
    AMP res;

    ca = comp_value(a);
    cr = comp_neg(ca);

    res = comp_lookup(cr);
    return res;
}

AMP
amp_add(AMP a, AMP b)
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

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_add(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_ADD, a, b, sylvan_false, res);
    }
    return res;
}

AMP
amp_sub(AMP a, AMP b)
{
    // special cases
    if (b == C_ZERO) return a;
    if (a == C_ZERO) return amp_neg(b);

    // check cache
    AMP res;
    if (CACHE_AMP_OPS) {
        if (cache_get3(CACHE_AMP_SUB, a, b, sylvan_false, &res)) {
            return res; // TODO: counters for these cache lookups/puts
        }
    }

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_sub(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_SUB, a, b, sylvan_false, res);
    }
    return res;
}

AMP
amp_mul(AMP a, AMP b)
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

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_mul(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_MUL, a, b, sylvan_false, res);
    }
    return res;
}

AMP
amp_div(AMP a, AMP b)
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
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_div(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put3(CACHE_AMP_DIV, a, b, sylvan_false, res);
    }
    return res;
}



/* Arithmetic operations on complex structs */

complex_t
comp_abs(complex_t a)
{
    complex_t res;
    res.r = sqrt( (a.r*a.r) + (a.i*a.i) );
    res.i = 0.0;
    return res;
}

complex_t
comp_cnj(complex_t a)
{
    complex_t res;
    res.r =  a.r;
    res.i = -a.i;
    return res;
}

complex_t
comp_neg(complex_t a)
{
    complex_t res;
    res.r = -a.r;
    res.i = -a.i;
    return res;
}

complex_t
comp_add(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

complex_t
comp_sub(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

complex_t
comp_mul(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

complex_t
comp_div(complex_t a, complex_t b)
{
    complex_t res;
    double denom;
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
comp_to_prob(complex_t a)
{
    double abs = sqrt( (a.r*a.r) + (a.i*a.i) );
    return (abs*abs);
}


/* Comparing complex values */

bool comp_exact_equal(complex_t a, complex_t b)
{
    return (a.r == b.r && a.i == b.i);
}

bool comp_approx_equal(complex_t a, complex_t b)
{
    return comp_epsilon_close(a, b, TOLERANCE);
}

bool comp_epsilon_close(complex_t a, complex_t b, double epsilon)
{
    return ( (fabs(a.r - b.r) < epsilon) && (fabs(a.i - b.i) < epsilon) );
}



/* Inserting / retrieving complex values from complex table */

AMP
comp_lookup(complex_t c)
{
    // TODO: catch comp_zero() / comp_one() here?
    uint64_t res;
    cmap_find_or_put(ctable, &c, &res);
    return (AMP) res;
}

complex_t
comp_value(AMP a)
{
    // special cases (lookup is read-only so this might make little/no difference)
    if (a == C_ZERO)    return comp_zero();
    if (a == C_ONE)     return comp_one();
    if (a == C_MIN_ONE) return comp_minus_one();

    // lookup
    complex_t * res;
    res = cmap_get(ctable, a);
    return *res;
}

/* used for gc of ctable */
complex_t
comp_value_old(AMP a)
{
    // special cases
    if (a == C_ZERO)    return comp_zero();
    if (a == C_ONE)     return comp_one();
    if (a == C_MIN_ONE) return comp_minus_one();

    complex_t * res;
    res = cmap_get(ctable_old, a);
    return *res;
}


/* Printing */

uint32_t comp_digits_print = 3;

void comp_print(complex_t c)
{
    comp_print_digits(c, comp_digits_print);
}

void comp_print_sci(complex_t c)
{
    comp_print_digits_sci(c, comp_digits_print);
}

void
comp_print_digits(complex_t c, uint32_t digits)
{
    if(c.r >= 0)
        printf(" ");
    printf("%.*f", digits, c.r);
    if (c.i > 0)
        printf("+%.*fi", digits, c.i);
    if (c.i < 0)
        printf("%.*fi", digits, c.i);
}

void
comp_print_digits_sci(complex_t c, uint32_t digits)
{
    if(c.r >= 0)
        printf(" ");
    printf("%.*e", digits, c.r);
    if (c.i > 0)
        printf("+%.*ei", digits, c.i);
    if (c.i < 0)
        printf("%.*ei", digits, c.i);
}

void
comp_print_bits(AMP a)
{
    print_bitvalues(ctable, a);
}



/* Managing the complex value table */

void
init_amplitude_table(size_t size)
{
    SIZE = size;
    ctable = cmap_create(SIZE);

    C_ONE     = comp_lookup(comp_one());
    C_ZERO    = comp_lookup(comp_zero());
    C_MIN_ONE = comp_lookup(comp_minus_one());

    Pi = 2.0 * acos(0.0);

    init_gates();
}

void
init_gates()
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
    gates[k][0] = C_ZERO; gates[k][1] = comp_lookup(comp_make(0.0, -1.0));
    gates[k][2] = comp_lookup(comp_make(0.0, 1.0));  gates[k][3] = C_ZERO;

    k = GATEID_Z;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = comp_lookup(comp_make(1.0/sqrt(2.0),0));
    gates[k][3] = comp_lookup(comp_make(-1.0/sqrt(2.0),0));

    k = GATEID_S;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(0.0, 1.0));

    k = GATEID_T;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/sqrt(2.0), 1.0/sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/sqrt(2.0), -1.0/sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5,-0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    k = GATEID_sqrtY;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(-0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5, 0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    init_phase_gates(255);
}

void
init_phase_gates(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    double angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / pow(2.0, (double) k);
        cartesian = comp_make_angle(angle);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = comp_lookup(cartesian);

        // backward rotation
        angle = -2*Pi / pow(2.0, (double) k);
        cartesian = comp_make_angle(angle);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = comp_lookup(cartesian);
    }
}

uint64_t
count_amplitude_table_enries()
{
    return cmap_count_entries(ctable);
}

void
free_amplitude_table()
{
    cmap_free(ctable);
}

void
init_new_empty_table()
{
    // point old to current (full) ctable
    ctable_old = ctable;

    // re-init new (empty) ctable
    init_amplitude_table(SIZE);
}

void
delete_old_table()
{
    // delete  old (full) table
    cmap_free(ctable_old);
}

AMP
move_from_old_to_new(AMP a)
{
    complex_t c = comp_value_old(a);
    return comp_lookup(c);
}
