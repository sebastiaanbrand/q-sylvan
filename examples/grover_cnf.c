#include "grover_cnf.h"

TASK_IMPL_5(QMDD, qmdd_grover_cnf_iteration, QMDD, qmdd, BDDVAR, n, BDDVAR, k, BDDVAR, clauses, int*, oracle)
{
    // Compute the results of each clause
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[clause*k+l] < 0) {
                qmdd = qmdd_gate(qmdd, GATEID_X, abs(oracle[clause*k+l])-1);
            }
        }
        // "oracle" call (Save results of clause in ancilla)
        qmdd = qmdd_cgate3(qmdd,GATEID_X,abs(oracle[clause*k])-1,abs(oracle[clause*k+1])-1,abs(oracle[clause*k+2])-1,n+1+clause);
        // Uncompute X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[clause*k+l] < 0) {
                qmdd = qmdd_gate(qmdd, GATEID_X, abs(oracle[clause*k+l])-1);
            }
        }
    }
    // Oracle flip (since the MCX is in the wrong direction, we use a Z and 2*H around the target)
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        qmdd = qmdd_gate(qmdd, GATEID_X, n+1+clause);
    }
    qmdd = qmdd_gate(qmdd, GATEID_H, n);
    qmdd = qmdd_cgate_range(qmdd,GATEID_Z,n,n+clauses-1,n+clauses);
    qmdd = qmdd_gate(qmdd, GATEID_H, n);
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        qmdd = qmdd_gate(qmdd, GATEID_X, n+1+clause);
    }
    // Uncompute all clauses
    for (BDDVAR clause = clauses; clause != 0; clause--) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[(clause-1)*k+l] < 0) {
                qmdd = qmdd_gate(qmdd, GATEID_X, abs(oracle[(clause-1)*k+l])-1);
            }
        }
        qmdd = qmdd_cgate3(qmdd,GATEID_X,abs(oracle[(clause-1)*k])-1,abs(oracle[(clause-1)*k+1])-1,abs(oracle[(clause-1)*k+2])-1,n+1+(clause-1));
        // Uncompute X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[(clause-1)*k+l] < 0) {
                qmdd = qmdd_gate(qmdd, GATEID_X, abs(oracle[(clause-1)*k+l])-1);
            }
        }
    }

    // H on all qubits (also ancilla to make CZ a CX)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qmdd = qmdd_gate(qmdd, GATEID_H, qubit);
    }

    // X on all qubits (except ancilla)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qmdd = qmdd_gate(qmdd, GATEID_X, qubit);
    }

    // Controlled-NOT over all qubits (on ancilla)
    qmdd = qmdd_cgate_range(qmdd, GATEID_X, 0, n-1, n);

    // X on all qubits (except ancilla)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qmdd = qmdd_gate(qmdd, GATEID_X, qubit);
    }

    // H on all qubits (also ancilla to maqubite CZ a CX)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qmdd = qmdd_gate(qmdd, GATEID_H, qubit);
    }

    return qmdd;
}

QMDD
qmdd_grover_cnf(BDDVAR n, int* oracle, BDDVAR k, BDDVAR clauses, BDDVAR n_answers)
{
    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) / n_answers ) );
    // start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(n+1+clauses) );
    for (BDDVAR qubit = 0; qubit < n+1+clauses; qubit++) init[qubit] = 0;
    // init[n] = 1; 
    QMDD qmdd = qmdd_create_basis_state(n+1+clauses, init);
    qmdd = qmdd_gate(qmdd, GATEID_X, n);

    // H on all qubits
    for (BDDVAR qubit = 0; qubit < n+1; qubit++) {
        qmdd = qmdd_gate(qmdd, GATEID_H, qubit);
    }

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        qmdd = qmdd_grover_cnf_iteration(qmdd, n, k, clauses, oracle);
    }

    free(init);

    return qmdd;
}
