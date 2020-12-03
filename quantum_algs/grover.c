#include "grover.h"
#include "sylvan_qdd_complex.h"

uint64_t
qdd_grover_approx_number_of_gates(BDDVAR nbits)
{
    uint64_t iterations = floor( 3.14159265359/4.0 * sqrt( pow(2,nbits) ) );
    uint64_t init_gates = (nbits + 1);
    uint64_t gts_per_it = (4 + 4*nbits);
    return (init_gates + iterations*gts_per_it);
}

bool *
qdd_grover_random_flag(BDDVAR nbits)
{
    bool *flag = malloc( (sizeof(bool)*nbits) );
    for (BDDVAR k = 0; k < nbits; k++) flag[k] = (rand() % 2);
    return flag;
}

static bool *
append_one(BDDVAR n, bool *flag)
{
    bool *res = malloc ( sizeof(bool)*(n+1) );
    for (BDDVAR k = 0; k < n; k++) res[k] = flag[k];
    res[n] = 1;
    return res;
}


TASK_IMPL_3(QDD, qdd_grover_iteration, QDD, qdd, BDDVAR, n, bool*, oracle)
{
    // "oracle" call  (apply -1 flag to desired amplitude)
    qdd = qdd_gate(qdd, GATEID_H, n);
    qdd = qdd_all_control_phase(qdd, n+1, oracle);
    qdd = qdd_gate(qdd, GATEID_H, n);

    // H on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_X, k);
    
    // Controlled Z over all qubits (except ancilla)
    qdd = qdd_cgate_range(qdd, GATEID_Z, 0, n-2, n-1);
    
    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_X, k);

    // H on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);
    return qdd;
}

QDD
qdd_grover(BDDVAR n, bool *flag)
{   
    LACE_ME;

    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(n+1) );
    for (BDDVAR k = 0; k < n; k++) init[k] = 0;
    init[n] = 1; 
    QDD qdd = qdd_create_basis_state(n+1, init);

    // H on all qubits
    for (BDDVAR k = 0; k < n+1; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    bool *oracle = append_one(n, flag);
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_iteration(qdd, n, oracle);
    }

    return qdd;
}


QDD
qdd_grover_matrix(BDDVAR n, bool *flag)
{
    LACE_ME;

    BDDVAR nqubits = n + 1;

    // Number of Grover iterations
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // Create matrix QDDs
    QDD all_H, first_n_H, first_n_X, oracle, first_n_CZ, grov_it;
    uint32_t *gate_list = malloc( sizeof(uint32_t)*(nqubits) );
    int *cgate_options = malloc( sizeof(int)*(nqubits) );

    // H on all qubits
    all_H = qdd_create_single_qubit_gates_same(nqubits, GATEID_H);

    // H on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_H;
    gate_list[n] = GATEID_I;
    first_n_H = qdd_create_single_qubit_gates(nqubits, gate_list);

    // oracle
    for (BDDVAR k = 0; k < n; k++) cgate_options[k] = flag[k];
    cgate_options[n] = 2; // ancilla qubit is target qubit of this multi-cgate
    oracle = qdd_create_multi_cgate(nqubits, cgate_options, GATEID_X);

    // X on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_X;
    gate_list[n] = GATEID_I;
    first_n_X = qdd_create_single_qubit_gates(nqubits, gate_list);

    // CZ gate on all qubits except ancilla
    for (BDDVAR k = 0; k < n-1; k++) cgate_options[k] = 1; // controls"
    cgate_options[n-1] = 2; // target last qubit before ancilla
    cgate_options[n] = -1; // do nothing to ancilla
    first_n_CZ = qdd_create_multi_cgate(nqubits, cgate_options, GATEID_Z);
    
    // Grover iteration = oracle + mean inversion
    grov_it = oracle;
    grov_it = qdd_matmat_mult(first_n_H,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_CZ, grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_H,  grov_it, nqubits);

    // Now, actually apply the circuit:
    // 1. Start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(nqubits) );
    for (BDDVAR k = 0; k < n; k++) init[k] = 0;
    init[n] = 1; 
    QDD state = qdd_create_basis_state(nqubits, init);

    // 2. H on all qubits
    state = qdd_matvec_mult(all_H, state, nqubits);

    // 3. Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        state = qdd_matvec_mult(grov_it, state, nqubits);
    }

    return state;
}
