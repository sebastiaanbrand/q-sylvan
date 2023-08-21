#ifndef WGT_COMPLEX_H
#define WGT_COMPLEX_H

#include "sylvan_edge_weights.h"
#include "edge_weight_storage/flt.h"


/******************<Implementation of edge_weights interface>******************/

complex_t *weight_complex_malloc();
void _weight_complex_value(void *wgt_store, AADD_WGT a, complex_t *res);
AADD_WGT weight_complex_lookup(complex_t *a);
AADD_WGT _weight_complex_lookup_ptr(complex_t *a, void *wgt_store);

void init_complex_one_zero(void *wgt_store);

void weight_complex_abs(complex_t *a);
void weight_complex_neg(complex_t *a);
void weight_complex_sqr(complex_t *a);
void weight_complex_add(complex_t *a, complex_t *b);
void weight_complex_sub(complex_t *a, complex_t *b);
void weight_complex_mul(complex_t *a, complex_t *b);
void weight_complex_div(complex_t *a, complex_t *b);
bool weight_complex_eq(complex_t *a, complex_t *b);
bool weight_complex_eps_close(complex_t *a, complex_t *b, double eps);
bool weight_complex_greater(complex_t *a, complex_t *b);

AADD_WGT wgt_complex_norm_L2(AADD_WGT *low, AADD_WGT *high);
AADD_WGT wgt_complex_get_low_L2normed(AADD_WGT high);

void weight_complex_fprint(FILE *stream, complex_t *a);


static inline AADD_WGT
complex_lookup_angle(fl_t theta, fl_t mag)
{
	complex_t c = cmake_angle(theta, mag);
	return weight_lookup(&c);
}

static inline AADD_WGT
complex_lookup(fl_t r, fl_t i)
{
	complex_t c = cmake(r, i);
	return weight_lookup(&c);
}

/*****************</Implementation of edge_weights interface>******************/

#endif
