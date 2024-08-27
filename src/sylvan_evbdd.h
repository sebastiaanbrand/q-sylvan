#ifndef SYLVAN_EVBDD_H
#define SYLVAN_EVBDD_H

#include <stdbool.h>
#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef uint64_t EVBDD;
typedef uint64_t EVBDD_WGT;  // Edge weight
typedef uint64_t EVBDD_TARG; // Edge target

static const EVBDD_TARG  EVBDD_TERMINAL = 1;
static const BDDVAR     EVBDD_INVALID_VAR = UINT8_MAX;

typedef enum weight_norm_strategy {
    NORM_LOW,
    NORM_MAX,
    NORM_MIN,
    NORM_L2,
    n_norm_strategies
} weight_norm_strategy_t;


/*******************<Garbage collection, references, marking>******************/

VOID_TASK_DECL_1(evbdd_gc_mark_rec, EVBDD);
#define evbdd_gc_mark_rec(a) RUN(evbdd_gc_mark_rec, a)

/**
 * Store the pointer <a> in the pointers table.
 */
void evbdd_protect(EVBDD* a);

/**
 * Delete the pointer <a> from the pointers table.
 */
void evbdd_unprotect(EVBDD* a);

/**
 * Compute the number of pointers in the pointers table.
 */
size_t evbdd_count_protected(void);

/**
 * Push an EVBDD variable to the pointer reference stack.
 * During gc the variable will be inspected and the contents will be marked.
 */
void evbdd_refs_pushptr(const EVBDD *a);

/**
 * Pop the last <amount> EVBDD variables from the pointer reference stack.
 */
void evbdd_refs_popptr(size_t amount);

/**
 * Push an EVBDD to the values reference stack.
 * During garbage collection the references EVBDD will be marked.
 */
EVBDD evbdd_refs_push(EVBDD a);

/**
 * Pop the last <amount> EVBDDs from the values reference stack.
 */
void evbdd_refs_pop(long amount);

/**
 * Push a Task that returns an EVBDD to the tasks reference stack.
 * Usage: evbdd_refs_spawn(SPAWN(function, ...));
 */
void evbdd_refs_spawn(Task *t);

/**
 * Pop a Task from the task reference stack.
 * Usage: EVBDD result = evbdd_refs_sync(SYNC(function));
 */
EVBDD evbdd_refs_sync(EVBDD a);

/******************</Garbage collection, references, marking>******************/





/*************************<Cleaning edge weight table>*************************/

/* enabled by default */
void evbdd_set_auto_gc_wgt_table(bool enabled);
/* default 0.5 */
void evbdd_set_gc_wgt_table_thres(double fraction_filled);
double evbdd_get_gc_wgt_table_thres();
void evbdd_gc_wgt_table();
bool evbdd_test_gc_wgt_table();

/**
 * Recursive function for moving weights from old to new edge weight table.
 */
#define _fill_new_wgt_table(a) (RUN(_fill_new_wgt_table, a))
TASK_DECL_1(EVBDD, _fill_new_wgt_table, EVBDD);

/************************</Cleaning edge weight table>*************************/





/******************************<Initialization>********************************/

/**
 * Similar initialization as for MTBDDs + edge weight table init.
 * Setting tolerance to -1 uses default tolerance.
 * real table: stores 2 real values per edge weight, instead of 1 tuple
 * NOTE: this function doesn't currently check if the combination of table
 * sizes (edge weight table + node table) works in combination with using 
 * a real-table or complex-table.
 * - init_wgt_tab_entries() is a pointer to a function which is called after gc 
 *   of edge table. This can be used to reinitialize edge weight table values 
 *   outside of any EVBDD. Can be NULL.
 */
void sylvan_init_evbdd(size_t min_wgt_tablesize, size_t max_wgt_tablesize, double wgt_tab_tolerance, int edge_weigth_backend, int norm_strat, void *init_wgt_tab_entries);
void sylvan_init_evbdd_defaults(size_t min_wgt_tablesize, size_t max_wgt_tablesize);
void evbdd_set_caching_granularity(int granularity);

/*****************************</Initialization>********************************/





/**************************<Matrix/vector operations>**************************/

/**
 * Get the EVBDD value corresponding to the given <path> (bool array of length
 * nvars(a)).
 */
EVBDD_WGT evbdd_getvalue(EVBDD a, bool* path);

/**
 * Recursive implementation of vector addition.
 */
#define evbdd_plus(a,b) (RUN(evbdd_plus,a,b))
TASK_DECL_2(EVBDD, evbdd_plus, EVBDD, EVBDD);


/* Computes Mat * |vec> (Wrapper function) */
#define evbdd_matvec_mult(mat,vec,nvars) (RUN(evbdd_matvec_mult,mat,vec,nvars))
TASK_DECL_3(EVBDD, evbdd_matvec_mult, EVBDD, EVBDD, BDDVAR);

/* Computes A*B (note generally AB != BA) (Wrapper function) */
#define evbdd_matmat_mult(a,b,nvars) (RUN(evbdd_matmat_mult,a,b,nvars))
TASK_DECL_3(EVBDD, evbdd_matmat_mult, EVBDD, EVBDD, BDDVAR);

/**
 * Recursive implementation of matrix-vector mult and matrix-matrix mult.
 */
TASK_DECL_4(EVBDD, evbdd_matvec_mult_rec, EVBDD, EVBDD, BDDVAR, BDDVAR);
TASK_DECL_4(EVBDD, evbdd_matmat_mult_rec, EVBDD, EVBDD, BDDVAR, BDDVAR);


/**
 * Computes inner product of two vectors <b|a> 
 * (Note that if b contains complex values, the complex conjugate is taken)
*/
#define evbdd_inner_product(a,b,nvars) (RUN(evbdd_inner_product,a,b,nvars,0))
TASK_DECL_4(EVBDD_WGT, evbdd_inner_product, EVBDD, EVBDD, BDDVAR, BDDVAR);

/**
 * Increases all the variable number in EVBDD a by k (used for tensor product)
 * 
 * @param a EVBDD over n vars {j, j+1, ..., j+n-1} (generally j=0)
 * 
 * @return EVBDD over n vars {j+k, j+k+1, ..., j+k+n-1}
 */
EVBDD evbdd_increase_all_vars(EVBDD a, int k);

/* Replace the terminal node in a with b (effectively stacks a and b) 
* (used for tensor product)
*/
EVBDD evbdd_replace_terminal(EVBDD a, EVBDD_TARG b);

/**
 * @param a EVBDD over vars 0...n-1 (n = nvars_a)
 * @param b EVBDD over vars 0...m-1
 * @param nvars_a number of vars of EVBDD a
 * 
 * @return EVBDD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
EVBDD evbdd_tensor_prod(EVBDD a, EVBDD b, BDDVAR nvars_a);

/**
 * Computes the tensor product of vec (tensor) vec
 * 
 * @param a EVBDD over vars 0...n-1 (n = nvars_a)
 * @param b EVBDD over vars 0...m-1
 * @param nvars_a number of vars of EVBDD a
 * 
 * @return EVBDD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
#define evbdd_vec_tensor_prod(a, b, nvars_a) evbdd_tensor_prod(a,b,nvars_a)

/**
 * Computes the tensor product of mat (tensor) mat
 * 
 * @param a EVBDD over vars 0...2n-1 (n = nvars_a)
 * @param b EVBDD over vars 0...2m-1
 * @param nvars_a number of vars of EVBDD a (counting only unprimed)
 * 
 * @return EVBDD over vars 0...2n-1...(2n+2m)-1, representing a (tensor) b
 */
#define evbdd_mat_tensor_prod(a, b, nvars_a) evbdd_tensor_prod(a,b,2*nvars_a)

/*************************</Matrix/vector operations>**************************/





/***************************<EVBDD utility functions>***************************/

/**
 * Count the number of EVBDD nodes.
 */
uint64_t evbdd_countnodes(EVBDD a);

/**************************</EVBDD utility functions>***************************/





/**************************<Printing & file writing>***************************/

/**
 * Write a .dot representation of a given EVBDD
 */
void evbdd_fprintdot(FILE *out, EVBDD a, bool draw_zeros);

/*************************</Printing & file writing>***************************/





/********************************<Debug stuff>*********************************/

/**
 * Checks if the vectors encoded by two EVBDDs are the same. In principle the
 * root edges of two identical EVBDDs should always be the same, so this is 
 * mostly a debug/testing function.
 * 
 * @param n Number of variables
 * @param exact If true, vec(a) should equal vec(b) exactly (float equal),
 * otherwise allow for preset float equivalence margin.
 * 
 * @returns True iff vec(a) == vec(b).
 */
bool evbdd_equivalent(EVBDD a, EVBDD b, int n, bool exact, bool verbose);

/** Sanity check to see if the EVBDD variables are ordered and < nvars. */
#define evbdd_is_ordered(a, nvars) evbdd_is_ordered_rec(a, 0, nvars)
bool evbdd_is_ordered_rec(EVBDD a, BDD nextvar, BDD nvars);

void evbdd_printnodes(EVBDD a);
bool _next_bitstring(bool *x, int n);
void _print_bitstring(bool *x, int n, bool backwards);
uint64_t bitarray_to_int(bool *x, int n, bool MSB_first);
bool * int_to_bitarray(uint64_t n, int length, bool MSB_first);
bool bit_from_int(uint64_t a, uint8_t index);
void reverse_bit_array(bool *x, int length);

/*******************************</Debug stuff>*********************************/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
