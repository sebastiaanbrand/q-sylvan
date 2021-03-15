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


/**********************************************
Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation
**********************************************/
//#define qdd_cos(fac,div) cosl((long double)(fac)*Pi/(long double)(div))
//#define qdd_sin(fac,div) sinl((long double)(fac)*Pi/(long double)(div))
static long double Pi;    // set value of global Pi

// Table parameters
static const double default_tolerance = 1e-14;
static double tolerance;
static amp_storage_backend_t amp_backend;
size_t table_size;

// Keep estimate for number of entries for gc purposes
size_t table_entries_est;
DECLARE_THREAD_LOCAL(table_entries_local, size_t); // these are added to _est
static const uint64_t table_entries_local_buffer = 1000; // every 1000 entries

// Actual table (old is used for gc purposes)
cmap_t *ctable;
cmap_t *ctable_old;

// TODO: "merge" ctable and rtable code in a cleaner way
//size_t rtable_entries_est;
//DECLARE_THREAD_LOCAL(rtable_entries_local, size_t); // these are added to _est
//static const uint64_t rtable_entries_local_buffer = 1000; // every 1000 entries
rmap_t *rtable;
rmap_t *rtable_old;

tree_map_t *rtree;
tree_map_t *rtree_old;

static bool CACHE_AMP_OPS = true;


/* Shorthand functions for making complex numbers */

complex_t
comp_make(fl_t r, fl_t i)
{
    complex_t res;
    res.r = r;
    res.i = i;
    return res;
}

complex_t
comp_make_angle(fl_t theta)
{
    complex_t c;
    c.r = flt_cos(theta);
    c.i = flt_sin(theta);
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

fl_t
comp_qmake(int a, int b, int c)
{
    return (( a +  b * flt_sqrt(2.0)) / c);
}


/* Cache puts / gets */

// deterministically order {a,b} for caching commuting op(a,b)
static void
order_inputs(AMP *a, AMP *b) 
{
    AMP x = (*a > *b) ? *a : *b;
    AMP y = (*a > *b) ? *b : *a;
    *a = x;
    *b = y;
}

static void
cache_put_add(AMP a, AMP b, AMP res)
{
    order_inputs(&a, &b);
    if (cache_put3(CACHE_AMP_ADD, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_ADD_CACHEDPUT);
    }
}

static bool
cache_get_add(AMP a, AMP b, AMP *res)
{
    order_inputs(&a, &b);
    if (cache_get3(CACHE_AMP_ADD, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_ADD_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_sub(AMP a, AMP b, AMP res)
{
    if (cache_put3(CACHE_AMP_SUB, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_SUB_CACHEDPUT);
    }
}

static bool
cache_get_sub(AMP a, AMP b, AMP *res)
{
    if (cache_get3(CACHE_AMP_SUB, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_SUB_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_mul(AMP a, AMP b, AMP res)
{
    order_inputs(&a, &b);
    if (cache_put3(CACHE_AMP_MUL, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_MUL_CACHEDPUT);
    }
    if (CACHE_INV_OPS) {
        // put inverse as well (empirically seems not so beneficial)
        if (cache_put3(CACHE_AMP_DIV, res, b, sylvan_false, a)) {
            sylvan_stats_count(AMP_DIV_CACHEDPUT);
        }
        if (cache_put3(CACHE_AMP_DIV, res, a, sylvan_false, b)) {
            sylvan_stats_count(AMP_DIV_CACHEDPUT);
        }
    }
}

static bool
cache_get_mul(AMP a, AMP b, AMP *res)
{
    order_inputs(&a, &b);
    if (cache_get3(CACHE_AMP_MUL, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_MUL_CACHED);
        return true;
    }
    return false;
}

static void
cache_put_div(AMP a, AMP b, AMP res)
{
    if (cache_put3(CACHE_AMP_DIV, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_DIV_CACHEDPUT);
    }
    if (CACHE_INV_OPS) {
        // put inverse as well (empirically seems beneficial)
        order_inputs(&b, &res);
        if (cache_put3(CACHE_AMP_MUL, b, res, sylvan_false, a)) {
            sylvan_stats_count(AMP_MUL_CACHEDPUT);
        }
    }
}

static bool
cache_get_div(AMP a, AMP b, AMP *res)
{
    if (cache_get3(CACHE_AMP_DIV, a, b, sylvan_false, res)) {
        sylvan_stats_count(AMP_DIV_CACHED);
        return true;
    }
    return false;
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
        if (cache_get_add(a, b, &res)) return res;
    }

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_add(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put_add(a, b, res);
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
        if (cache_get_sub(a, b, &res)) return res;
    }

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_sub(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put_sub(a, b, res);
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
        if (cache_get_mul(a, b, &res)) return res;
    }

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_mul(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put_mul(a, b, res);
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
        if (cache_get_div(a, b, &res)) return res;
    }

    // compute and hash result to ctable
    complex_t ca, cb, cr;
    ca = comp_value(a);
    cb = comp_value(b);
    cr = comp_div(ca, cb);
    res = comp_lookup(cr);

    // insert in cache
    if (CACHE_AMP_OPS) {
        cache_put_div(a, b, res);
    }
    return res;
}

double
amp_to_prob(AMP a)
{
    return comp_to_prob(comp_value(a));
}

AMP
prob_to_amp(double a)
{
    complex_t c;
    c.r = flt_sqrt(a);
    c.i = 0;
    return comp_lookup(c);
}


/* Arithmetic operations on complex structs */

complex_t
comp_abs(complex_t a)
{
    complex_t res;
    res.r = flt_sqrt( (a.r*a.r) + (a.i*a.i) );
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
    fl_t denom;
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
    double abs = flt_sqrt ( (a.r*a.r) + (a.i*a.i) );
    return (abs*abs);
}


/* Comparing complex values */

bool comp_exact_equal(complex_t a, complex_t b)
{
    return (a.r == b.r && a.i == b.i);
}

bool comp_approx_equal(complex_t a, complex_t b)
{
    return comp_epsilon_close(a, b, cmap_get_tolerance());
}

bool comp_epsilon_close(complex_t a, complex_t b, double epsilon)
{
    return ( (flt_abs(a.r - b.r) < epsilon) && (flt_abs(a.i - b.i) < epsilon) );
}

/* Comparing AMPs */

bool amp_exact_equal(AMP a, AMP b)
{
    return comp_exact_equal(comp_value(a), comp_value(b));
}

bool amp_approx_equal(AMP a, AMP b)
{
    return comp_approx_equal(comp_value(a), comp_value(b));
}

bool amp_epsilon_close(AMP a, AMP b, double epsilon)
{
    return comp_epsilon_close(comp_value(a), comp_value(b), epsilon);
}


/* normalization of two amps */

AMP
amp_normalize_low(AMP *low, AMP *high)
{
    // Normalize using low if low != 0
    AMP norm;
    if(*low != C_ZERO){
        *high = amp_div(*high, *low);
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
amp_normalize_largest(AMP *low, AMP *high)
{
    AMP norm;
    if (*low == *high) {
        norm  = *low;
        *low  = C_ONE;
        *high = C_ONE;
        return norm;
    }

    // Normalize using the absolute greatest value
    complex_t cl = comp_value(*low);
    complex_t ch = comp_value(*high);
    if ( (cl.r*cl.r + cl.i*cl.i)  >=  (ch.r*ch.r + ch.i*ch.i) ) {
        *high = amp_div(*high, *low);
        norm = *low;
        *low  = C_ONE;
    }
    else {
        *low = amp_div(*low, *high);
        norm  = *high;
        *high = C_ONE;
    }
    return norm;
}




/* Inserting / retrieving complex values from complex table */

AMP
comp_lookup(complex_t c)
{
    // TODO: catch comp_zero() / comp_one() here?
    uint64_t res;
    bool success;
    res = comp_try_lookup(c, &success);
    if (!success) {
        printf("Amplitude table full!\n");
        exit(1);
    }
    return (AMP) res;
}

static AMP
comp_try_lookup_ctable(complex_t c, bool *success)
{
    uint64_t res;
    int present = cmap_find_or_put(ctable, &c, &res);
    if (present == 0) {
        *success = true;
        table_entries_local += 1;
        if (table_entries_local >= table_entries_local_buffer) {
            __sync_fetch_and_add(&table_entries_est, table_entries_local);
            table_entries_local = 0;
        }
    }
    else if (present == 1) *success = true;
    else                   *success = false;
    return (AMP) res;
}

static AMP
try_lookup_rtable(double a, bool *success) // TODO: change to fl_t
{
    uint64_t res;
    int present = rmap_find_or_put(rtable, &a, &res);
    if (present == 0) {
        *success = true;
        table_entries_local += 1;
        if (table_entries_local >= table_entries_local_buffer) {
            __sync_fetch_and_add(&table_entries_est, table_entries_local);
            table_entries_local = 0;
        }
    }
    else if (present == 1) *success = true;
    else                   *success = false;
    return (AMP) res;
}

static AMP
try_lookup_rtree(fl_t a, bool *success)
{
    unsigned int res;
    int present = tree_map_find_or_put(rtree, a, &res);
    if (present == 0) {
        *success = true;
        table_entries_local += 1;
        if (table_entries_local >= table_entries_local_buffer) {
            __sync_fetch_and_add(&table_entries_est, table_entries_local);
            table_entries_local = 0;
        }
    }
    else if (present == 1) *success = true;
    else                   *success = false;
    return (AMP) res;
}

static AMP
pack_indices_rtable(AMP r, AMP i)
{
    uint64_t index_size = (int) ceil(log2(table_size));
    AMP res = ((r << index_size) | i);
    return res;
}

static void
unpack_indices_rtable(AMP bundle, AMP *index_r, AMP *index_i)
{
    uint64_t index_size = (int) ceil(log2(table_size));
    *index_r = (bundle >> index_size);
    *index_i = (bundle & ((1<<index_size)-1));
}

static AMP
comp_try_lookup_rtable(complex_t c, bool *success)
{
    bool suc1, suc2;
    AMP amp1 = try_lookup_rtable(c.r, &suc1);
    AMP amp2 = try_lookup_rtable(c.i, &suc2);
    *success = (suc1 && suc2);
    return pack_indices_rtable(amp1, amp2);
}

static AMP
comp_try_lookup_rtree(complex_t c, bool *success)
{
    bool suc1, suc2;
    AMP amp1 = try_lookup_rtree(c.r, &suc1);
    AMP amp2 = try_lookup_rtree(c.i, &suc2);
    *success = (suc1 && suc2);
    return pack_indices_rtable(amp1, amp2);
}

AMP
comp_try_lookup(complex_t c, bool *success)
{
    if (amp_backend == REAL_HASHMAP) {
        return comp_try_lookup_rtable(c, success);
    }
    else if (amp_backend == COMP_HASHMAP) {
        return comp_try_lookup_ctable(c, success);
    }
    else if (amp_backend == REAL_TREE) {
        return comp_try_lookup_rtree(c, success);
    }
    else {
        printf("lookup: backend not recognized\n");
        exit(1);
    }
}

complex_t
comp_value(AMP a)
{
    // special cases (lookup is read-only so this might make little/no difference)
    if (a == C_ZERO)    return comp_zero();
    if (a == C_ONE)     return comp_one();
    if (a == C_MIN_ONE) return comp_minus_one();

    // lookup
    complex_t res;
    if (amp_backend == REAL_HASHMAP) {
        AMP idx_r, idx_i;
        double *r, *i;
        unpack_indices_rtable(a, &idx_r, &idx_i);
        r = rmap_get(rtable, idx_r);
        i = rmap_get(rtable, idx_i);
        res.r = *r;
        res.i = *i;
    }
    else if (amp_backend == COMP_HASHMAP) {
        res = *cmap_get(ctable, a);
    }
    else if (amp_backend == REAL_TREE) {
        AMP idx_r, idx_i;
        fl_t *r, *i;
        unpack_indices_rtable(a, &idx_r, &idx_i);
        r = tree_map_get(rtree, idx_r);
        i = tree_map_get(rtree, idx_i);
        res.r = *r;
        res.i = *i;
    }
    else {
        printf("get value: backend not recognized\n");
        exit(1);
    }
    return res;
}

/* used for gc of ctable */
complex_t
comp_value_old(AMP a)
{
    // special cases
    if (a == C_ZERO)    return comp_zero();
    if (a == C_ONE)     return comp_one();
    if (a == C_MIN_ONE) return comp_minus_one();

    complex_t res;
    if (amp_backend == REAL_HASHMAP) {
        AMP idx_r, idx_i;
        double *r, *i;
        unpack_indices_rtable(a, &idx_r, &idx_i);
        r = rmap_get(rtable_old, idx_r);
        i = rmap_get(rtable_old, idx_i);
        res.r = *r;
        res.i = *i;
    }
    else if (amp_backend == COMP_HASHMAP) {
        res = *cmap_get(ctable_old, a);
    }
    else if (amp_backend == REAL_TREE) {
        AMP idx_r, idx_i;
        fl_t *r, *i;
        unpack_indices_rtable(a, &idx_r, &idx_i);
        r = tree_map_get(rtree_old, idx_r);
        i = tree_map_get(rtree_old, idx_i);
        res.r = *r;
        res.i = *i;
    }
    else {
        printf("get value old: backend not recognized\n");
        exit(1);
    }
    return res;
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
    printf("%.*Lf", digits, (long double) c.r);
    if (c.i > 0)
        printf("+%.*Lfi", digits, (long double) c.i);
    if (c.i < 0)
        printf("%.*Lfi", digits, (long double) c.i);
}

void
comp_print_digits_sci(complex_t c, uint32_t digits)
{
    if(c.r >= 0)
        printf(" ");
    printf("%.*Le", digits, (long double) c.r);
    if (c.i > 0)
        printf("+%.*Lei", digits, (long double) c.i);
    if (c.i < 0)
        printf("%.*Lei", digits, (long double) c.i);
}

void
comp_print_bits(AMP a)
{
    print_bitvalues(ctable, a);
}


/*************************** <dynamic custom gates> ***************************/

uint32_t next_custom_id; // set to 0 in init

uint32_t
get_custom_gate_id()
{
    next_custom_id++;
    if (next_custom_id >= num_dynamic_gates) {
        // max custom gates used, reset ID counter to 0 and clear opcache
        next_custom_id = 0;
        LACE_ME;
        sylvan_clear_cache();
    }
    return num_static_gates + next_custom_id; // index offset by num_static_gates
}

uint32_t
GATEID_Rz(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    double theta_over_2 = Pi * a;
    AMP u00, u11;
    u00 = comp_lookup(comp_make_angle(-theta_over_2));
    u11 = comp_lookup(comp_make_angle(theta_over_2));
    gates[gate_id][0] = u00;    gates[gate_id][1] = C_ZERO;
    gates[gate_id][2] = C_ZERO; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Rx(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    fl_t theta_over_2 = Pi * a;
    AMP u00, u01, u10, u11;
    u00 = comp_lookup(comp_make(flt_cos(theta_over_2), 0.0));
    u01 = comp_lookup(comp_make(0.0, -flt_sin(theta_over_2)));
    u10 = comp_lookup(comp_make(0.0, -flt_sin(theta_over_2)));
    u11 = comp_lookup(comp_make(flt_cos(theta_over_2), 0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Ry(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    fl_t theta_over_2 = Pi * a;
    AMP u00, u01, u10, u11;
    u00 = comp_lookup(comp_make(flt_cos(theta_over_2),  0.0));
    u01 = comp_lookup(comp_make(-flt_sin(theta_over_2), 0.0));
    u10 = comp_lookup(comp_make(flt_sin(theta_over_2),  0.0));
    u11 = comp_lookup(comp_make(flt_cos(theta_over_2),  0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

/************************** </dynamic custom gates> ***************************/


/* Managing the complex value table */

void
init_amplitude_table(size_t size, long double tol, amp_storage_backend_t backend)
{
    tolerance = (tol < 0) ? default_tolerance : tol;
    table_size = size;
    amp_backend = backend;

    table_entries_est = 0;
    table_entries_local = 0;

    // TODO: something cleaner than this if structure
    if (amp_backend == REAL_HASHMAP) {
        rtable = rmap_create(table_size, tolerance);
    }
    else if (amp_backend == COMP_HASHMAP) {
        ctable = cmap_create(table_size, tolerance);
    }
    else if (amp_backend == REAL_TREE) {
        rtree = tree_map_create(table_size, tolerance);
        //printf("real tree not integrated into complex value handling yet\n");
        //exit(1);
    }
    else {
        printf("init amp table: backend not recognized\n");
        exit(1);
    }
    
    // NOTE: the sum of the local counters sometimes exceeds the actual total
    // number of entries (when just counting the global value with atomic adds
    // after every insert). This might be because 'ctable_entries_local = 0' 
    // doesn't set it to 0 for all threads. Since it's only an estimate it is
    // not a huge issue though.
    // TODO: figure out how to get lace to handle this counting better


    C_ONE     = comp_lookup(comp_one());
    C_ZERO    = comp_lookup(comp_zero());
    C_MIN_ONE = comp_lookup(comp_minus_one());

    Pi = 2.0 * flt_acos(0.0);

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
    gates[k][0] = gates[k][1] = gates[k][2] = comp_lookup(comp_make(1.0/flt_sqrt(2.0),0));
    gates[k][3] = comp_lookup(comp_make(-1.0/flt_sqrt(2.0),0));

    k = GATEID_S;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(0.0, 1.0));

    k = GATEID_Sdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(0.0, -1.0));

    k = GATEID_T;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5,-0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    k = GATEID_sqrtY;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(-0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5, 0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    init_phase_gates(255);

    next_custom_id = 0;
}

void
init_phase_gates(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    fl_t angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / (fl_t)(1<<k);
        cartesian = comp_make_angle(angle);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = comp_lookup(cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
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

uint64_t
get_table_entries_estimate()
{
    return table_entries_est;
}

uint64_t
get_table_size()
{
    return table_size;
}

void
free_amplitude_table()
{
    if (amp_backend == REAL_HASHMAP) {
        rmap_free(rtable);
    }
    else if (amp_backend == COMP_HASHMAP) {
        cmap_free(ctable);
    }
    else if (amp_backend == REAL_TREE) {
        tree_map_free(rtree);
    }
}

void
init_new_empty_table()
{
    // point old to current (full) ctable
    if (amp_backend == REAL_HASHMAP) {
        rtable_old = rtable;
    }
    else if (amp_backend == COMP_HASHMAP) {
        ctable_old = ctable;
    }
    else if (amp_backend == REAL_TREE) {
        rtree_old = rtree;
    }

    // re-init new (empty) ctable
    double tolerance = 0;
    if (amp_backend == REAL_HASHMAP) {
        tolerance = rmap_get_tolerance();
    }
    else if (amp_backend == COMP_HASHMAP) {
        tolerance = cmap_get_tolerance();
    }
    else if (amp_backend == REAL_TREE) {
        tolerance = tree_map_get_tolerance();
    }
    init_amplitude_table(table_size, tolerance, amp_backend);
}

void
delete_old_table()
{
    // delete  old (full) table
    if (amp_backend == REAL_HASHMAP) {
        rmap_free(rtable_old);
    }
    else if (amp_backend == COMP_HASHMAP) {
        cmap_free(ctable_old);
    }
    else if (amp_backend == REAL_TREE) {
        tree_map_free(rtree_old);
    }
}

AMP
move_from_old_to_new(AMP a)
{
    complex_t c = comp_value_old(a);
    return comp_lookup(c);
}
