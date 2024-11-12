#include <stdio.h>
#include <stdlib.h>
#include "sylvan_edge_weights_complex.h"


/**********************<Some static utility functions>*************************/

static void
comp_cart_to_polar(fl_t r, fl_t i, fl_t *magnitude, fl_t *angle)
{
    *magnitude = flt_sqrt(r*r + i*i);
    if (r != 0) {
        *angle = flt_atan2(i, r);
    }
    else {
        if (i > 0) *angle = flt_acos(0.0); //  pi/2
        else if (i < 0) *angle = -flt_acos(0.0); // -pi/2
        else *angle = 0; // i = 0
    }
}

/**********************</Some static utility functions>************************/




/*****************<Implementation of edge_weights interface>*******************/

complex_t *
weight_complex_malloc()
{
    complex_t *res = malloc(sizeof(complex_t));
    return res;
}

void
_weight_complex_value(void *wgt_store, EVBDD_WGT a, complex_t *res)
{
    if (a == EVBDD_ZERO)         *res = czero();
    else if (a == EVBDD_ONE)     *res = cone();
    else if (a == EVBDD_MIN_ONE) *res = cmone();
    *res = *(complex_t*)(wgt_store_get(wgt_store, a)); // ?
}

EVBDD_WGT
_weight_complex_lookup_ptr(complex_t *a, void *wgt_store)
{
    // TODO: catch czero() / cone() here?
    uint64_t res;
    bool success;

    int present = wgt_store_find_or_put(wgt_store, a, &res);
    if (present == -1) {
        success = false;
    } else if (present == 0) { 
        success = true;
        wgt_table_gc_inc_entries_estimate();
    } else {
        success = true;
    }

    if (!success) {
        fprintf(stderr, "Amplitude table full!\n");
        exit(1);
    }
    return (EVBDD_WGT) res; 
}

EVBDD_WGT
weight_complex_lookup(complex_t *a)
{
    return _weight_complex_lookup_ptr(a, wgt_storage);
}

void
init_complex_one_zero(void *wgt_store)
{
    complex_t a;
    a = cone();     EVBDD_ONE     = _weight_complex_lookup_ptr(&a, wgt_store);
    a = czero();    EVBDD_ZERO    = _weight_complex_lookup_ptr(&a, wgt_store);
    a = cmone();    EVBDD_MIN_ONE = _weight_complex_lookup_ptr(&a, wgt_store);
}

void
weight_complex_abs(complex_t *a)
{
    a->r = flt_sqrt( (a->r*a->r) + (a->i*a->i) );
    a->i = 0.0;
}

void
weight_complex_neg(complex_t *a)
{
    a->r = -(a->r);
    a->i = -(a->i);
}

void
weight_complex_conj(complex_t *a)
{
    a->i = -(a->i);
}

void
weight_complex_sqr(complex_t *a)
{
    weight_complex_mul(a, a);
}

void
weight_complex_add(complex_t *a, complex_t *b)
{
    a->r = a->r + b->r;
    a->i = a->i + b->i;
}

void
weight_complex_sub(complex_t *a, complex_t *b)
{
    a->r = a->r - b->r;
    a->i = a->i - b->i;
}

void
weight_complex_mul(complex_t *a, complex_t *b)
{
    complex_t tmp;
    tmp.r = a->r * b->r - a->i * b->i;
    tmp.i = a->r * b->i + a->i * b->r;
    a->r = tmp.r;
    a->i = tmp.i;
}

void
weight_complex_div(complex_t *a, complex_t *b)
{
    complex_t tmp;
    fl_t denom;
    if (b->i == 0.0) {
        tmp.r = a->r / b->r;
        tmp.i = a->i / b->r;
    } else {
        denom = b->r * b->r + b->i * b->i;
        tmp.r = (a->r * b->r + a->i * b->i) / denom;
        tmp.i = (a->i * b->r - a->r * b->i) / denom;
    }
    a->r = tmp.r;
    a->i = tmp.i;
}

bool
weight_complex_eq(complex_t *a, complex_t *b)
{
    return ( (a->r == b->r) && (a->i == b->i) );
}

bool
weight_complex_eps_close(complex_t *a, complex_t *b, double eps)
{
    return ( (flt_abs(a->r - b->r) <= eps) && (flt_abs(a->i - b->i) <= eps) );
}

bool
weight_complex_greater(complex_t *a, complex_t *b)
{
    return ( (a->r*a->r + a->i*a->i) > (b->r*b->r + b->i*b->i) );
}

EVBDD_WGT
wgt_complex_norm_L2(EVBDD_WGT *low, EVBDD_WGT *high)
{
    // normalize such that |low|^2 + |high|^2 = 1, and low in R+

    // Deal with cases where one weight is 0 (both 0 shouldn't end up here)
    if (*low == EVBDD_ZERO) {
        EVBDD_WGT res = *high;
        *high = EVBDD_ONE;
        return res;
    }
    else if (*high == EVBDD_ZERO){
        EVBDD_WGT res = *low;
        *low = EVBDD_ONE;
        return res;
    }

    // TODO: add caching ?

    complex_t a, b;
    weight_value(*low, &a);
    weight_value(*high, &b);

    // convert to polar form
    fl_t mag_a, mag_b, theta_a, theta_b;
    comp_cart_to_polar(a.r, a.i, &mag_a, &theta_a);
    comp_cart_to_polar(b.r, b.i, &mag_b, &theta_b);

    // normalize magnitudes
    fl_t _norm = flt_sqrt(mag_a*mag_a + mag_b*mag_b);
    mag_a = mag_a / _norm;
    mag_b = mag_b / _norm;

    // normalize phase (subtract theta_a from both) to have low \in R+
    theta_b = theta_b - theta_a;
    // theta_a will be set to 0 for a', but needed for norm

    // convert to cartesian form
    a = cmake(mag_a, 0); // theta_a = 0
    b = cmake_angle(theta_b, mag_b);
    complex_t c_norm = cmake_angle(theta_a, _norm);

    // return
    *low  = weight_complex_lookup(&a);
    *high = weight_complex_lookup(&b);
    return weight_complex_lookup(&c_norm);
}

EVBDD_WGT
wgt_complex_get_low_L2normed(EVBDD_WGT high)
{
    // Get low from high, assuming |low|^2 + |high|^2 = 1, and low \in R+:
    // a = sqrt(1 - |b|^2)
    if (high == EVBDD_ZERO) return EVBDD_ONE;
    if (high == EVBDD_ONE || high == EVBDD_MIN_ONE) return EVBDD_ZERO;

    complex_t b;
    weight_value(high, &b);
    fl_t mag_b = (b.r*b.r + b.i*b.i);
    if (mag_b > 1.0) {
        if (mag_b > 1.0 + 1e-6 && mag_b <= 1.1) { 
            printf("Warning: |b| = %.15lf > 1.0\n", mag_b);
            printf("Continuing with |b| = 1.0\n");
        } else if (mag_b > 1.1) {
            printf("Value error in L2 norm: |b| = %.15lf > 1.0\n", mag_b);
            exit(1);
        }
        mag_b = 1.0;
    }
    fl_t a = flt_sqrt(1.0 - mag_b);

    complex_t c = cmake(a, 0);
    return weight_complex_lookup(&c);
}

void
weight_complex_fprint(FILE *stream, complex_t *a)
{
    int digits = 3;
    if(a->r >= 0)
        fprintf(stream, " ");
    fprintf(stream, "%.*Lf", digits, (long double) a->r);
    if (a->i > 0)
        fprintf(stream, "+%.*Lfi", digits, (long double) a->i);
    if (a->i < 0)
        fprintf(stream, "%.*Lfi", digits, (long double) a->i);
}

/*****************</Implementation of edge_weights interface>******************/
