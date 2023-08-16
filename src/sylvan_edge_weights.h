#ifndef SYLVAN_EDGE_WEIGHTS_H
#define SYLVAN_EDGE_WEIGHTS_H

#include <stdint.h>
#include <stdio.h>
#include <edge_weight_storage/wgt_storage_interface.h>

typedef uint64_t AADD_WGT;  // AADD edge weights (indices to table entries)

extern AADD_WGT AADD_ONE;
extern AADD_WGT AADD_ZERO;
extern AADD_WGT AADD_MIN_ONE;

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

extern void sylvan_init_edge_weights(size_t size, double tol, edge_weight_type_t edge_weight_type, wgt_storage_backend_t backend);
extern void init_edge_weight_functions(edge_weight_type_t edge_weight_type);
extern void init_edge_weight_storage(size_t size, double tol, wgt_storage_backend_t backend, void **wgt_store);
extern void (*init_wgt_table_entries)(); // set by sylvan_init_aadd

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
extern AADD_WGT wgt_table_gc_keep(AADD_WGT a);

/************************</GC of edge weight table>****************************/





/******************<Interface for different edge_weight_types>*****************/

typedef weight_t (*weight_malloc_f)();
typedef void (*_weight_value_f)(void *wgt_store, AADD_WGT a, weight_t res);
typedef AADD_WGT (*weight_lookup_f)(weight_t a);
typedef AADD_WGT (*_weight_lookup_ptr_f)(weight_t a, void *wgt_store);

typedef void (*init_one_zero_f)(void *wgt_store);

/* Arithmetic operations on edge weights */
typedef void (*weight_abs_f)(weight_t a); // a <-- |a|
typedef void (*weight_neg_f)(weight_t a); // a <-- -a
typedef void (*weight_sqr_f)(weight_t a); // a <-- a^2
typedef void (*weight_add_f)(weight_t a, weight_t b); // a <-- a + b
typedef void (*weight_sub_f)(weight_t a, weight_t b); // a <-- a - b
typedef void (*weight_mul_f)(weight_t a, weight_t b); // a <-- a * b
typedef void (*weight_div_f)(weight_t a, weight_t b); // a <-- a / b
typedef bool (*weight_eq_f)(weight_t a, weight_t b); // returns true iff a == b
typedef bool (*weight_eps_close_f)(weight_t a, weight_t b, double eps); // returns true iff dist(a,b) < eps
typedef bool (*weight_greater_f)(weight_t a, weight_t b); // returns true iff a > b

/* Normalization methods */
typedef AADD_WGT (*wgt_norm_L2_f)(AADD_WGT *low, AADD_WGT *high);
typedef AADD_WGT (*wgt_get_low_L2normed_f)(AADD_WGT high);

typedef void (*weight_fprint_f)(FILE *stream, weight_t a);



extern weight_malloc_f 		weight_malloc;
extern _weight_value_f 		_weight_value;

extern weight_lookup_f 		weight_lookup;
extern _weight_lookup_ptr_f	_weight_lookup_ptr;
extern init_one_zero_f 		init_one_zero;
extern weight_abs_f 		weight_abs;
extern weight_neg_f 		weight_neg;
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





/*********************<Arithmetic functions on AADD_WGT's>*********************/

/* Arithmetic operations on AADD_WGT's */
AADD_WGT wgt_abs(AADD_WGT a); // returns |a|
AADD_WGT wgt_neg(AADD_WGT a); // returns -a
AADD_WGT wgt_add(AADD_WGT a, AADD_WGT b); // returns a + b
AADD_WGT wgt_sub(AADD_WGT a, AADD_WGT b); // returns a - b
AADD_WGT wgt_mul(AADD_WGT a, AADD_WGT b); // returns a * b
AADD_WGT wgt_div(AADD_WGT a, AADD_WGT b); // returns a / b

/********************</Arithmetic functions on AADD_WGT's>*********************/





/*************************<Comparators on AADD_WGT's>**************************/

bool wgt_eq(AADD_WGT a, AADD_WGT b);
bool wgt_eps_close(AADD_WGT a, AADD_WGT b, double eps);
bool wgt_approx_eq(AADD_WGT a, AADD_WGT b);

/************************</Comparators on AADD_WGT's>**************************/





/*************************<Edge weight normalization>**************************/

AADD_WGT wgt_norm_low(AADD_WGT *low, AADD_WGT *high);
AADD_WGT wgt_norm_largest(AADD_WGT *low, AADD_WGT *high);
// wgt_norm_L2() is in the interface because it's too complicated to implement
// without assumptions on the underlying data type of the edge weights.

/************************</Edge weight normalization>**************************/





/************************<Printing & utility functions>************************/

void wgt_fprint(FILE *stream, AADD_WGT a);

/************************<Printing & utility functions>************************/

#endif // SYLVAN_EDGE_WEIGHTS_H
