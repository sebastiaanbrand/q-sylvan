#include <stdbool.h>

#include "sylvan.h"

/* Random bit array of lenght 'nqubits' */
// bool *qdd_grover_cnf_random_flag(BDDVAR nqubits, BDDVAR clauses);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QDD.
 */
QDD qdd_grover_cnf(BDDVAR n, bool* oracle, BDDVAR clauses, BDDVAR n_answers);
#define qdd_grover_cnf_iteration(qdd,n,clauses,oracle) (CALL(qdd_grover_cnf_iteration,qdd,n,clauses,oracle));
TASK_DECL_4(QDD, qdd_grover_cnf_iteration, QDD, BDDVAR, BDDVAR, bool*);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
// QDD qdd_grover_cnf_matrix(BDDVAR n, bool *flag, BDDVAR n_answers);
