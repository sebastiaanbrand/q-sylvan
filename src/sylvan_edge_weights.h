#ifndef SYLVAN_EDGE_WEIGHTS_H
#define SYLVAN_EDGE_WEIGHTS_H

#include <stdint.h>
#include <stdio.h>
#include <edge_weight_storage/wgt_storage_interface.h>

typedef uint64_t EVBDD_WGT;  // EVBDD edge weights (indices to table entries)

extern EVBDD_WGT EVBDD_ONE;
extern EVBDD_WGT EVBDD_ZERO;
extern EVBDD_WGT EVBDD_MIN_ONE;

typedef void *weight_t;

typedef enum edge_weight_type {
    WGT_DOUBLE,
    WGT_COMPLEX_128,
    WGT_RATIONAL_128,
    n_wgt_types
} edge_weight_type_t;


// Actual table (new used for gc purposes)
extern void *wgt_storage; // TODO: move to source file?
extern void *wgt_storage_new;


/**********************<Managing the edge weight table>************************/

extern void sylvan_init_edge_weights(size_t min_tablesize, size_t max_tablesize, double tol, edge_weight_type_t edge_weight_type, wgt_storage_backend_t backend);
extern void init_edge_weight_functions(edge_weight_type_t edge_weight_type);
extern void init_edge_weight_storage(size_t size, double tol, wgt_storage_backend_t backend, void **wgt_store);
extern void (*init_wgt_table_entries)(); // set by sylvan_init_evbdd

extern uint64_t sylvan_get_edge_weight_table_size();
extern double sylvan_edge_weights_tolerance();
extern uint64_t sylvan_edge_weights_count_entries();
extern void sylvan_edge_weights_free();

/*********************</Managing the edge weight table>************************/





/*************************<GC of edge weight table>****************************/

extern void init_edge_weight_storage_gc();
extern uint64_t wgt_table_entries_estimate();
extern void wgt_table_gc_inc_entries_estimate();
extern void wgt_table_gc_init_new(void (*init_wgt_table_entries)());
extern void wgt_table_gc_delete_old();
extern EVBDD_WGT wgt_table_gc_keep(EVBDD_WGT a);

/************************</GC of edge weight table>****************************/





/******************<Interface for different edge_weight_types>*****************/

typedef weight_t (*weight_malloc_f)();
typedef void (*_weight_value_f)(void *wgt_store, EVBDD_WGT a, weight_t res);
typedef EVBDD_WGT (*weight_lookup_f)(weight_t a);
typedef EVBDD_WGT (*_weight_lookup_ptr_f)(weight_t a, void *wgt_store);

typedef void (*init_one_zero_f)(void *wgt_store);

/* Arithmetic operations on edge weights */
typedef void (*weight_abs_f)(weight_t a); // a <-- |a|
typedef void (*weight_neg_f)(weight_t a); // a <-- -a
typedef void (*weight_conj_f)(weight_t a); // a <-- a*
typedef void (*weight_sqr_f)(weight_t a); // a <-- a^2
typedef void (*weight_add_f)(weight_t a, weight_t b); // a <-- a + b
typedef void (*weight_sub_f)(weight_t a, weight_t b); // a <-- a - b
typedef void (*weight_mul_f)(weight_t a, weight_t b); // a <-- a * b
typedef void (*weight_div_f)(weight_t a, weight_t b); // a <-- a / b
typedef bool (*weight_eq_f)(weight_t a, weight_t b); // returns true iff a == b
typedef bool (*weight_eps_close_f)(weight_t a, weight_t b, double eps); // returns true iff dist(a,b) < eps
typedef bool (*weight_greater_f)(weight_t a, weight_t b); // returns true iff |a| > |b|

/* Normalization methods */
typedef EVBDD_WGT (*wgt_norm_L2_f)(EVBDD_WGT *low, EVBDD_WGT *high);
typedef EVBDD_WGT (*wgt_get_low_L2normed_f)(EVBDD_WGT high);

typedef void (*weight_fprint_f)(FILE *stream, weight_t a);



extern weight_malloc_f 		weight_malloc;
extern _weight_value_f 		_weight_value;

extern weight_lookup_f 		weight_lookup;
extern _weight_lookup_ptr_f	_weight_lookup_ptr;
extern init_one_zero_f 		init_one_zero;
extern weight_abs_f 		weight_abs;
extern weight_neg_f 		weight_neg;
extern weight_conj_f        weight_conj;
extern weight_sqr_f 		weight_sqr;
extern weight_add_f 		weight_add;
extern weight_sub_f 		weight_sub;
extern weight_mul_f 		weight_mul;
extern weight_div_f 		weight_div;
extern weight_eq_f 			weight_eq;
extern weight_eps_close_f 	weight_eps_close;
extern weight_greater_f		weight_greater;

extern wgt_norm_L2_f		wgt_norm_L2;
extern wgt_get_low_L2normed_f		wgt_get_low_L2normed;

extern weight_fprint_f 		weight_fprint;


#define weight_lookup_ptr(a) _weight_lookup_ptr(a, wgt_storage)
#define weight_value(a, res) _weight_value(wgt_storage, a, res)
#define weight_approx_eq(a,b) weight_eps_close(a,b,sylvan_edge_weights_tolerance())

/*****************</Interface for different edge_weight_types>*****************/





/********************<For caching arithmetic operations>***********************/

void wgt_set_inverse_chaching(bool on);

/*******************</For caching arithmetic operations>***********************/





/*********************<Arithmetic functions on EVBDD_WGT's>*********************/

/* Arithmetic operations on EVBDD_WGT's */
EVBDD_WGT wgt_abs(EVBDD_WGT a); // returns |a|
EVBDD_WGT wgt_neg(EVBDD_WGT a); // returns -a
EVBDD_WGT wgt_conj(EVBDD_WGT a); // returns a*
EVBDD_WGT wgt_add(EVBDD_WGT a, EVBDD_WGT b); // returns a + b
EVBDD_WGT wgt_sub(EVBDD_WGT a, EVBDD_WGT b); // returns a - b
EVBDD_WGT wgt_mul(EVBDD_WGT a, EVBDD_WGT b); // returns a * b
EVBDD_WGT wgt_div(EVBDD_WGT a, EVBDD_WGT b); // returns a / b

/********************</Arithmetic functions on EVBDD_WGT's>*********************/





/*************************<Comparators on EVBDD_WGT's>**************************/

bool wgt_eq(EVBDD_WGT a, EVBDD_WGT b);
bool wgt_eps_close(EVBDD_WGT a, EVBDD_WGT b, double eps);
bool wgt_approx_eq(EVBDD_WGT a, EVBDD_WGT b);

/************************</Comparators on EVBDD_WGT's>**************************/





/*************************<Edge weight normalization>**************************/

EVBDD_WGT wgt_norm_low(EVBDD_WGT *low, EVBDD_WGT *high);
EVBDD_WGT wgt_norm_max(EVBDD_WGT *low, EVBDD_WGT *high);
EVBDD_WGT wgt_norm_min(EVBDD_WGT *low, EVBDD_WGT *high);
// wgt_norm_L2() is in the interface because it's too complicated to implement
// without assumptions on the underlying data type of the edge weights.

/************************</Edge weight normalization>**************************/





/************************<Printing & utility functions>************************/

void wgt_fprint(FILE *stream, EVBDD_WGT a);

/************************<Printing & utility functions>************************/

#endif // SYLVAN_EDGE_WEIGHTS_H
