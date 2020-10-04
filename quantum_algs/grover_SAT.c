#include "grover.h"
#include "sylvan_qdd_complex.h"

bool *
qdd_grover_cnf_random_flag(BDDVAR nqubits, BDDVAR clauses)
{
    bool *flag = malloc( (sizeof(bool)*nqubits*clauses) );
    for (BDDVAR k = 0; k < nqubits*clauses; k++)
        flag[k] = 0;//(rand() % 2);
    return flag;
}


TASK_IMPL_3(QDD, qdd_grover_cnf_iteration, QDD, qdd, BDDVAR, n, bool*, flag, BDDVAR, clauses)
{

    bool *clause = malloc( (sizeof(bool)*n) );
    // each clause is a flag on its own, loop over all flags
    for (BDDVAR k = 0; k < clauses; k++) {
        for (BDDVAR l = 0; l < n; l++)
            clause[l] = flag[k*n+l];
        // "oracle" call  (apply -1 flag to desired amplitude)
        qdd = qdd_all_control_phase(qdd, n, clause);
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
    QDD all_H, clause_oracle, oracle, phase, grov_it;
    all_H  = qdd_create_single_qubit_gates_same(n, GATEID_H);
    phase  = qdd_grover_cnf_phase_flips(n);

    // Oracle matrix QDD is a little bit more difficult
    // Create for each clause a matrix QDD and multiply with the previous matrix QDD
    // The endresult will be one matrix QDD with all clause oracles in it
    bool *clause = malloc( (sizeof(bool)*n) );
    for (BDDVAR l = 0; l < n; l++)
        clause[l] = flag[l];
    oracle = qdd_create_all_control_phase(n, clause);
    // each clause is a flag on its own, loop over all flags
    for (BDDVAR k = 1; k < clauses; k++) {
        for (BDDVAR l = 0; l < n; l++)
            clause[l] = flag[k*n+l];
        // "oracle" call  (apply -1 flag to desired amplitude)
        clause_oracle = qdd_create_all_control_phase(n, clause);
        oracle = qdd_matmat_mult(clause_oracle, oracle, n);
    }

    // Grover iteration matrix G = H * Phase * H * oracle
    // (oracle on the left because we apply G|\psi> from right to left)
    grov_it = oracle;
    grov_it = qdd_matmat_mult(all_H, grov_it, n);
    grov_it = qdd_matmat_mult(phase, grov_it, n);
    grov_it = qdd_matmat_mult(all_H, grov_it, n);

    // Actual circuit, start with all 0 state
    QDD state = qdd_create_all_zero_state(n);

    // H on all qubits
    state = qdd_matvec_mult(all_H, state, n);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        state = qdd_matvec_mult(grov_it, state, n);
    }

    return state;
}
