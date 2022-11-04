#include <qsylvan.h>

/* default = false */
void qmdd_grover_set_verbose(bool v);

/* Approximate number of gates for 'nbits'-bit Grover */
uint64_t qmdd_grover_approx_number_of_gates(BDDVAR nbits);

/* Random bit array of lenght 'nbits'. IMPORTANT: this returns a malloc array. */
bool *qmdd_grover_random_flag(BDDVAR nbits);

/* Bit array 11..1 of lenght 'nbits'. IMPORTANT: this returns a malloc array. */
bool *qmdd_grover_ones_flag(BDDVAR nbits);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QMDD.
 */
QMDD qmdd_grover(BDDVAR n, bool* flag);
#define qmdd_grover_iteration(qmdd,n,oracle) (CALL(qmdd_grover_iteration,qmdd,n,oracle))
TASK_DECL_3(QMDD, qmdd_grover_iteration, QMDD, BDDVAR, bool*);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QMDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
QMDD qmdd_grover_matrix(BDDVAR n, bool *flag);

/**
 * Instead of applying iteration matrix G, R times, applies G^t, R/t times.
 * For t=1 this function acts as the "normal" matrix implementation of Grover.
 * (matrix is returned as arg for nodecount purposes)
 */
QMDD qmdd_grover_matrix_multi_its(BDDVAR n, bool *flag, int t, QMDD *matrix);
