#include <stdbool.h>

#include "sylvan.h"

/**
 * Exectutes Grover on an n qubit circuit with a single flagged element.
 * 
 * @param n Number of qubits.
 * @param flag Binary representation of some integer in [0, 2^n]
 */
QDD qdd_grover(BDDVAR n, bool* flag);
#define qdd_grover_iteration(qdd,n,flag) (CALL(qdd_grover_iteration,qdd,n,flag));
TASK_DECL_3(QDD, qdd_grover_iteration, QDD, BDDVAR, bool*);
