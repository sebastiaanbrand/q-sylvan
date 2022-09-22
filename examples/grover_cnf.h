#include <qsylvan.h>

/* Random bit array of lenght 'nqubits' */
// bool *qmdd_grover_cnf_random_flag(BDDVAR nqubits, BDDVAR clauses);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QMDD.
 */
QMDD qmdd_grover_cnf(BDDVAR n, int* oracle, BDDVAR k, BDDVAR clauses, BDDVAR n_answers);
#define qmdd_grover_cnf_iteration(qmdd,n,k,clauses,oracle) (CALL(qmdd_grover_cnf_iteration,qmdd,n,k,clauses,oracle));
TASK_DECL_5(QMDD, qmdd_grover_cnf_iteration, QMDD, BDDVAR, BDDVAR,BDDVAR, int*);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QMDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
// QMDD qmdd_grover_cnf_matrix(BDDVAR n, bool *flag, BDDVAR n_answers);
