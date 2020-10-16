#include "grover_cnf.h"
#include "sylvan_qdd_complex.h"

TASK_IMPL_4(QDD, qdd_grover_cnf_iteration, QDD, qdd, BDDVAR, n, BDDVAR, clauses, bool*, oracle)
{
    // Compute the results of each clause
    for (BDDVAR k = 0; k < clauses; k++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < n; l++) {
            if(!oracle[k*n+l]) qdd = qdd_gate(qdd, GATEID_X, l);
        }
        // "oracle" call (Save results of clause in ancilla)
        qdd = qdd_cgate_range(qdd,GATEID_X,0,n-1,n+1+k);
        // Uncompute X gates
        for (BDDVAR l = 0; l < n; l++) {
            if(!oracle[k*n+l]) qdd = qdd_gate(qdd, GATEID_X, l);
        }
    }
    // Oracle flip (since the MCX is in the wrong direction, we use a Z and 2*H around the target)
    for (BDDVAR k = 0; k < clauses; k++) qdd = qdd_gate(qdd, GATEID_X, n+1+k);
    qdd = qdd_gate(qdd, GATEID_H, n);
    qdd = qdd_cgate_range(qdd,GATEID_Z,n,n+clauses-1,n+clauses);
    qdd = qdd_gate(qdd, GATEID_H, n);
    for (BDDVAR k = 0; k < clauses; k++) qdd = qdd_gate(qdd, GATEID_X, n+1+k);
    // Uncompute all clauses
    for (BDDVAR k = 0; k < clauses; k++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < n; l++) {
            if(!oracle[(clauses-1-k)*n+l]) qdd = qdd_gate(qdd, GATEID_X, l);
        }
        qdd = qdd_cgate_range(qdd,GATEID_X,0,n-1,n+1+(clauses-1-k));
        // Uncompute X gates
        for (BDDVAR l = 0; l < n; l++)
            if(!oracle[(clauses-1-k)*n+l]) qdd = qdd_gate(qdd, GATEID_X, l);
    }


    // H on all qubits (also ancilla to make CZ a CX)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_X, k);

    // Controlled-NOT over all qubits (on ancilla)
    qdd = qdd_cgate_range(qdd, GATEID_X, 0, n-1, n);

    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_X, k);

    // H on all qubits (also ancilla to make CZ a CX)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    return qdd;
}

QDD
qdd_grover_cnf(BDDVAR n, bool* oracle, BDDVAR clauses, BDDVAR n_answers)
{
    LACE_ME;
    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) / n_answers ) );
    // start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(n+1+clauses) );
    for (BDDVAR k = 0; k < n+1+clauses; k++) init[k] = 0;
    // init[n] = 1; 

    QDD qdd = qdd_create_basis_state(n+1+clauses, init);
    qdd = qdd_gate(qdd, GATEID_X, n);

    // H on all qubits
    for (BDDVAR k = 0; k < n+1; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_cnf_iteration(qdd, n, clauses, oracle);
    }

    return qdd;
}
