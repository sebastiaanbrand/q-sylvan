#ifndef SYLVAN_EDGE_WEIGHTS_H
#define SYLVAN_EDGE_WEIGHTS_H

#include <stdint.h>
#include <stdio.h>
#include <edge_weight_storage/wgt_storage_interface.h>

typedef uint64_t AADD_WGT;  // AADD edge weights (indices to table entries)

extern AADD_WGT AADD_ONE;
extern AADD_WGT AADD_ZERO;
extern AADD_WGT AADD_MIN_ONE;

typedef void* weight_t;

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

extern weight_t (*weight_malloc)();
extern void (*_weight_value)(); // weight_value(WGT a, weight_type * res)
extern AADD_WGT (*weight_lookup)(); // weight_lookup_(weight_type a)
extern AADD_WGT (*_weight_lookup_ptr)(); // _weight_lookup_ptr(weight_type *a, void *wgt_store)
#define weight_lookup_ptr(a) _weight_lookup_ptr(a, wgt_storage)
#define weight_value(a, res) _weight_value(wgt_storage, a, res)

extern void (*init_one_zero)();

/* Arithmetic operations on edge weights */
extern void (*weight_abs)(); // a <-- |a|
extern void (*weight_neg)(); // a <-- -a
extern void (*weight_sqr)(); // a <-- a^2
extern void (*weight_add)(); // a <-- a + b
extern void (*weight_sub)(); // a <-- a - b
extern void (*weight_mul)(); // a <-- a * b
extern void (*weight_div)(); // a <-- a / b
extern bool (*weight_eq)(); // returns true iff a == b
extern bool (*weight_eps_close)(); // returns true iff dist(a,b) < eps
#define weight_approx_eq(a,b) weight_eps_close(a,b,sylvan_edge_weights_tolerance())
extern bool (*weight_greater)(); // returns true iff a > b

/* Normalization methods */
extern AADD_WGT (*wgt_norm_L2)(); // wgt_norm_L2(AADD_WGT *low, AADD_WGT *high)
extern AADD_WGT (*wgt_get_low_L2normed)(); // wgt_get_low_L2normed(AADD_WGT high)

extern void (*weight_fprint)(); // weight_fprint(FILE *stream, weight_t a)

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
