#ifndef SYLVAN_AADD_H
#define SYLVAN_AADD_H

#include <stdbool.h>
#include <stdint.h>


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef uint64_t AADD;
typedef uint64_t AADD_WGT;  // Edge weight
typedef uint64_t AADD_TARG; // Edge target

static const AADD_TARG  AADD_TERMINAL = 1;
static const BDDVAR     AADD_INVALID_VAR = UINT8_MAX;

typedef enum weight_norm_strategy {
    NORM_LOW,
    NORM_MAX,
    NORM_MIN,
    NORM_L2,
    n_norm_strategies
} weight_norm_strategy_t;


/*******************<Garbage collection, references, marking>******************/

VOID_TASK_DECL_1(aadd_gc_mark_rec, AADD);
#define aadd_gc_mark_rec(a) RUN(aadd_gc_mark_rec, a)

/**
 * Store the pointer <a> in the pointers table.
 */
void aadd_protect(AADD* a);

/**
 * Delete the pointer <a> from the pointers table.
 */
void aadd_unprotect(AADD* a);

/**
 * Compute the number of pointers in the pointers table.
 */
size_t aadd_count_protected(void);

/**
 * Push an AADD variable to the pointer reference stack.
 * During gc the variable will be inspected and the contents will be marked.
 */
void aadd_refs_pushptr(const AADD *a);

/**
 * Pop the last <amount> AADD variables from the pointer reference stack.
 */
void aadd_refs_popptr(size_t amount);

/**
 * Push an AADD to the values reference stack.
 * During garbage collection the references AADD will be marked.
 */
AADD aadd_refs_push(AADD a);

/**
 * Pop the last <amount> AADDs from the values reference stack.
 */
void aadd_refs_pop(long amount);

/**
 * Push a Task that returns an AADD to the tasks reference stack.
 * Usage: aadd_refs_spawn(SPAWN(function, ...));
 */
void aadd_refs_spawn(Task *t);

/**
 * Pop a Task from the task reference stack.
 * Usage: AADD result = aadd_refs_sync(SYNC(function));
 */
AADD aadd_refs_sync(AADD a);

/******************</Garbage collection, references, marking>******************/





/*************************<Cleaning edge weight table>*************************/

/* enabled by default */
void aadd_set_auto_gc_wgt_table(bool enabled);
/* default 0.5 */
void aadd_set_gc_wgt_table_thres(double fraction_filled);
double aadd_get_gc_wgt_table_thres();
void aadd_gc_wgt_table();
bool aadd_test_gc_wgt_table();

/**
 * Recursive function for moving weights from old to new edge weight table.
 */
#define _fill_new_wgt_table(a) (RUN(_fill_new_wgt_table, a))
TASK_DECL_1(AADD, _fill_new_wgt_table, AADD);

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
 *   outside of any AADD. Can be NULL.
 */
void sylvan_init_aadd(size_t min_wgt_tablesize, size_t max_wgt_tablesize, double wgt_tab_tolerance, int edge_weigth_backend, int norm_strat, void *init_wgt_tab_entries);
void sylvan_init_aadd_defaults(size_t min_wgt_tablesize, size_t max_wgt_tablesize);
void aadd_set_caching_granularity(int granularity);

/*****************************</Initialization>********************************/





/**************************<Matrix/vector operations>**************************/

/**
 * Get the AADD value corresponding to the given <path> (bool array of length
 * nvars(a)).
 */
AADD_WGT aadd_getvalue(AADD a, bool* path);

/**
 * Recursive implementation of vector addition.
 */
#define aadd_plus(a,b) (RUN(aadd_plus,a,b))
TASK_DECL_2(AADD, aadd_plus, AADD, AADD);


/* Computes Mat * |vec> (Wrapper function) */
#define aadd_matvec_mult(mat,vec,nvars) (RUN(aadd_matvec_mult,mat,vec,nvars))
TASK_DECL_3(AADD, aadd_matvec_mult, AADD, AADD, BDDVAR);

/* Computes A*B (note generally AB != BA) (Wrapper function) */
#define aadd_matmat_mult(a,b,nvars) (RUN(aadd_matmat_mult,a,b,nvars))
TASK_DECL_3(AADD, aadd_matmat_mult, AADD, AADD, BDDVAR);

/**
 * Recursive implementation of matrix-vector mult and matrix-matrix mult.
 */
TASK_DECL_4(AADD, aadd_matvec_mult_rec, AADD, AADD, BDDVAR, BDDVAR);
TASK_DECL_4(AADD, aadd_matmat_mult_rec, AADD, AADD, BDDVAR, BDDVAR);


/**
 * Computes inner product of two vectors <b|a> 
 * (Note that if b contains complex values, the complex conjugate is taken)
*/
#define aadd_inner_product(a,b,nvars) (RUN(aadd_inner_product,a,b,nvars,0))
TASK_DECL_4(AADD_WGT, aadd_inner_product, AADD, AADD, BDDVAR, BDDVAR);

/**
 * Increases all the variable number in AADD a by k (used for tensor product)
 * 
 * @param a AADD over n vars {j, j+1, ..., j+n-1} (generally j=0)
 * 
 * @return AADD over n vars {j+k, j+k+1, ..., j+k+n-1}
 */
AADD aadd_increase_all_vars(AADD a, int k);

/* Replace the terminal node in a with b (effectively stacks a and b) 
* (used for tensor product)
*/
AADD aadd_replace_terminal(AADD a, AADD_TARG b);

/**
 * @param a AADD over vars 0...n-1 (n = nvars_a)
 * @param b AADD over vars 0...m-1
 * @param nvars_a number of vars of AADD a
 * 
 * @return AADD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
AADD aadd_tensor_prod(AADD a, AADD b, BDDVAR nvars_a);

/**
 * Computes the tensor product of vec (tensor) vec
 * 
 * @param a AADD over vars 0...n-1 (n = nvars_a)
 * @param b AADD over vars 0...m-1
 * @param nvars_a number of vars of AADD a
 * 
 * @return AADD over vars 0...n-1...(n+m)-1, representing a (tensor) b
 */
#define aadd_vec_tensor_prod(a, b, nvars_a) aadd_tensor_prod(a,b,nvars_a)

/**
 * Computes the tensor product of mat (tensor) mat
 * 
 * @param a AADD over vars 0...2n-1 (n = nvars_a)
 * @param b AADD over vars 0...2m-1
 * @param nvars_a number of vars of AADD a (counting only unprimed)
 * 
 * @return AADD over vars 0...2n-1...(2n+2m)-1, representing a (tensor) b
 */
#define aadd_mat_tensor_prod(a, b, nvars_a) aadd_tensor_prod(a,b,2*nvars_a)

/*************************</Matrix/vector operations>**************************/





/***************************<AADD utility functions>***************************/

/**
 * Count the number of AADD nodes.
 */
uint64_t aadd_countnodes(AADD a);

/**************************</AADD utility functions>***************************/





/**************************<Printing & file writing>***************************/

/**
 * Write a .dot representation of a given AADD
 */
void aadd_fprintdot(FILE *out, AADD a, bool draw_zeros);

/*************************</Printing & file writing>***************************/





/********************************<Debug stuff>*********************************/

/**
 * Checks if the vectors encoded by two AADDs are the same. In principle the
 * root edges of two identical AADDs should always be the same, so this is 
 * mostly a debug/testing function.
 * 
 * @param n Number of variables
 * @param exact If true, vec(a) should equal vec(b) exactly (float equal),
 * otherwise allow for preset float equivalence margin.
 * 
 * @returns True iff vec(a) == vec(b).
 */
bool aadd_equivalent(AADD a, AADD b, int n, bool exact, bool verbose);

/** Sanity check to see if the AADD variables are ordered and < nvars. */
#define aadd_is_ordered(a, nvars) aadd_is_ordered_rec(a, 0, nvars)
bool aadd_is_ordered_rec(AADD a, BDD nextvar, BDD nvars);

void aadd_printnodes(AADD a);
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
