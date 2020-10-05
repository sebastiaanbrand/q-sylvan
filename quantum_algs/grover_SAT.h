#include <stdbool.h>

#include "sylvan.h"

/* Random bit array of lenght 'nqubits' */
bool *qdd_grover_cnf_random_flag(BDDVAR nqubits);

/**
 * Implementation of Grover where the gates are seen as functions applied to
 * the QDD.
 */
QDD qdd_grover_cnf(BDDVAR n, bool* flag, BDDVAR clauses, BDDVAR n_answers);
#define qdd_grover_cnf_iteration(qdd,n,flag,clauses) (CALL(qdd_grover_cnf_iteration,qdd,n,flag,clauses));
TASK_DECL_3(QDD, qdd_grover_cnf_iteration, QDD, BDDVAR, bool*, BDDVAR);

/**
 * Implementation of Grover where both the state vector and the gates are
 * represented as QDDs, and matrix-vector / matrix-matrix multiplication is
 * used to compute the result.
 */
QDD qdd_grover_cnf_matrix(BDDVAR n, BDDVAR clauses, bool *flag, BDDVAR n_answers);
