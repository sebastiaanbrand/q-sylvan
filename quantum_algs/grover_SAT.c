#include "grover.h"
#include "sylvan_qdd_complex.h"

bool *
qdd_grover_cnf_random_flag(BDDVAR nqubits, BDDVAR nclauses)
{
    // #solutions = pow(2, pow(2, nqubits) - nclauses - 1)
    bool *flag = malloc( (sizeof(bool)*nqubits*nclauses) );
    // Generate predefined k-SAT (generating random k-SAT where #answers is known is too much work)
    BDDVAR clauses = nclauses;
    // If #clauses > #options for n qubits, no solutions remain
    // In this case, copy already used clauses to satisfy nclauses
    if (nclauses >= pow(2, nqubits)) {
        clauses = pow(2, nqubits) - 1;
        for (BDDVAR k = clauses; k < nclauses; k++) {
            for (BDDVAR l = 0; l < nqubits; l++) {
                if (k >= pow(2,l))
                    flag[k*nqubits+l] = 0;
    }
    for (BDDVAR k = 0; k < clauses; k++) {
        BDDVAR m = k;
        for (BDDVAR l = 0; l < nqubits; l++) {
            flag[k*nqubits+l] = m & 1;
            m /= 2;
        }
    }
    return flag;
}

TASK_IMPL_3(QDD, qdd_grover_cnf_iteration, QDD, qdd, BDDVAR, n, bool*, flag, BDDVAR, clauses)
{
    // Compute the results of each clause
    for (BDDVAR k = 0; k < clauses; k++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < n; l++)
            if(!flag[k*n+l]) { qdd = qdd_gate(qdd, GATEID_X, l); }
        // "oracle" call (Save results of clause in ancilla)
        qdd = qdd_gate(qdd, GATEID_H, n+k);
        qdd = qdd_cgate_range(qdd,GATEID_Z,0,n,n+k);
        qdd = qdd_gate(qdd, GATEID_H, n+k);
        // Uncompute X gates
        for (BDDVAR l = 0; l < n; l++)
            if(!flag[k*n+l]) { qdd = qdd_gate(qdd, GATEID_X, l); }
    }
    for (BDDVAR k = 0; k < clauses; k++)
        qdd = qdd_gate(qdd, GATEID_X, n+k);
    qdd = qdd_cgate_range(qdd,GATEID_Z,0,n+k-1,n+k);
    for (BDDVAR k = 0; k < clauses; k++)
        qdd = qdd_gate(qdd, GATEID_X, n+k);
    // Uncompute all clauses
    for (BDDVAR k = 0; k < clauses; k++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < n; l++)
            if(!flag[k*n+l]) { qdd = qdd_gate(qdd, GATEID_X, l); }
        // "oracle" call (Save results of clause in ancilla)
        qdd = qdd_gate(qdd, GATEID_H, n+k);
        qdd = qdd_cgate_range(qdd,GATEID_Z,0,n,n+k);
        qdd = qdd_gate(qdd, GATEID_H, n+k);
        // Uncompute X gates
        for (BDDVAR l = 0; l < n; l++)
            if(!flag[k*n+l]) { qdd = qdd_gate(qdd, GATEID_X, l); }
    }

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++)
        qdd = qdd_gate(qdd, GATEID_H, k);

    // Phase on all amplitudes except |000...0>
    bool x[n];
    for(BDDVAR k = 0; k < n; k++)
        x[k] = 0;
    qdd = qdd_all_control_phase(qdd, n, x);
    qdd = qdd_scalar_mult(qdd, comp_minus_one());

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++)
        qdd = qdd_gate(qdd, GATEID_H, k);

    return qdd;
}

QDD
qdd_grover_cnf(BDDVAR n, bool* flag, BDDVAR clauses, BDDVAR n_answers)
{   
    LACE_ME;

    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) / n_answers ) );

    // start with all zero state |000...0>
    QDD qdd = qdd_create_all_zero_state(n);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_cnf_iteration(qdd, n, flag, clauses);
    }

    return qdd;
}

QDD
qdd_grover_cnf_phase_flips(BDDVAR n)
{
    // flip all the phases except for |00..0>
    bool x[n]; 
    for(BDDVAR k = 0; k < n; k++) x[k] = 0;
    QDD phase = qdd_create_all_control_phase(n, x);
    phase = qdd_scalar_mult(phase, comp_minus_one());
    return phase;
}

QDD
qdd_grover_cnf_matrix(BDDVAR n, BDDVAR clauses, bool *flag, BDDVAR n_answers)
{
    LACE_ME;

    // Number of Grover iterations
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) / n_answers ) );

    // Create matrix QDDs
    QDD all_H, all_I, oracle, oracle_half, oracle_flip, phase, grov_it;
    all_I  = qdd_create_single_qubit_gates_same(n+clauses, GATEID_I);
    all_H  = qdd_create_single_qubit_gates_same(n, GATEID_H);
    phase  = qdd_grover_cnf_phase_flips(n);





    // for (BDDVAR k = 0; k < clauses; k++)
    //     qdd = qdd_gate(qdd, GATEID_X, n+k);
    // qdd = qdd_cgate_range(qdd,GATEID_Z,0,n+k-1,n+k);
    // for (BDDVAR k = 0; k < clauses; k++)
    //     qdd = qdd_gate(qdd, GATEID_X, n+k);
    // // Uncompute all clauses

    // Oracle matrix QDD is a little bit more difficult
    // Create for each clause a matrix QDD and multiply with the previous matrix QDD
    // The endresult will be one matrix QDD with all clause oracles in it
    oracle_half = all_I;
    for (BDDVAR k = 0; k < clauses; k++) {
        // Make a matrix which contains all gates that can be put in parallel
        // e.g. take clause 3, which is (x1 V x2)
        // Transform (x1 V x2) to ~(~x1 ^ ~x2)
        QDD clause_oracle, ancilla_store;
        // Put a Pauli X gate on corresponding qubits (e.g. x1 and x2) and I on the others
        for (BDDVAR l = 0; l < n; l++)
            if(!flag[k*n+l]) { clause_oracle = qdd_gate(clause_oracle, GATEID_X, l); }
        // Save results of clause in ancilla (e.g. ancilla 3)
        // Place on this qubit a Hadamard (since we use control-Z) and I on the others
        for (BDDVAR l = 0; l < clauses; l++)
            if (l == k) { clause_oracle = qdd_gate(clause_oracle, GATEID_H, n+l); }

        // Prepare clause constraints
        oracle_half = qdd_matmat_mult(clause_oracle, oracle_half, n+clauses);
        // Store in ancilla
        ancilla_store = qdd_cgate_range(ancilla_store,GATEID_Z,0,n,n+k);
        oracle_half = qdd_matmat_mult(ancilla_store, oracle_half, n+clauses)
        // Undo clause constraints (uncompute)
        oracle_half = qdd_matmat_mult(clause_oracle, oracle_half, n+clauses);
    }
    // Oracle flip
    for (BDDVAR k = 0; k < clauses; k++)
        oracle_flip = qdd_gate(oracle_flip, GATEID_X, n+k);
    oracle_flip = qdd_cgate_range(oracle_flip,GATEID_Z,0,n+k-1,n+k);
    for (BDDVAR k = 0; k < clauses; k++)
        oracle_flip = qdd_gate(oracle_flip, GATEID_X, n+k);

    // Construct Oracle
    oracle = oracle_half;
    oracle = qdd_matmat_mult(oracle_flip, oracle, n+clauses);
    oracle = qdd_matmat_mult(oracle_half, oracle, n+clauses);

    // Grover iteration matrix G = H * Phase * H * oracle
    // (oracle on the left because we apply G|\psi> from right to left)
    grov_it = oracle;
    grov_it = qdd_matmat_mult(all_H, grov_it, n+clauses);
    grov_it = qdd_matmat_mult(phase, grov_it, n+clauses);
    grov_it = qdd_matmat_mult(all_H, grov_it, n+clauses);

    // Actual circuit, start with all 0 state
    QDD state = qdd_create_all_zero_state(n+clauses);

    // H on all qubits
    state = qdd_matvec_mult(all_H, state, n);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        state = qdd_matvec_mult(grov_it, state, n+clauses);
    }

    return state;
}
