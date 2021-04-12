#include <stdbool.h>

#include "sylvan.h"

/* default = false */
void qdd_grover_set_verbose(bool v);

/* Approximate number of gates for 'nbits'-bit Grover */
uint64_t qdd_grover_approx_number_of_gates(BDDVAR nbits);

/* Random bit array of lenght 'nbits' */
bool *qdd_grover_random_flag(BDDVAR nbits);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QDD.
 */
QDD qdd_grover(BDDVAR n, bool* flag);
#define qdd_grover_iteration(qdd,n,oracle) (CALL(qdd_grover_iteration,qdd,n,oracle))
TASK_DECL_3(QDD, qdd_grover_iteration, QDD, BDDVAR, bool*);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
QDD qdd_grover_matrix(BDDVAR n, bool *flag);

/**
 * Instead of applying iteration matrix G, R times, applies G^t, R/t times.
 * For t=1 this function acts as the "normal" matrix implementation of Grover.
 */
QDD qdd_grover_matrix_multi_its(BDDVAR n, bool *flag, int t);
