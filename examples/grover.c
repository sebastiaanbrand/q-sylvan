#include "grover.h"

static bool VERBOSE = false;
static bool EXP_BY_SQUARING = true;

void
qdd_grover_set_verbose(bool v)
{
    VERBOSE = v;
}

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
    free(init);

    // H on all qubits
    for (BDDVAR k = 0; k < n+1; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    bool *oracle = append_one(n, flag);
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_iteration(qdd, n, oracle);
    }
    free(oracle);

    return qdd;
}

QDD
qdd_grover_matrix(BDDVAR n, bool *flag)
{
    QDD m;
    return qdd_grover_matrix_multi_its(n, flag, 1, &m);
}

QDD qdd_grover_matrix_multi_its(BDDVAR n, bool *flag, int t, QDD *matrix)
{
    LACE_ME;

    BDDVAR nqubits = n + 1;

    // Number of Grover iterations
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // Create matrix QDDs
    QDD all_H, first_n_H, first_n_X, oracle, first_n_CZ, grov_it, grov_its, state;
    uint32_t *gate_list = malloc( sizeof(uint32_t)*(nqubits) );
    int *cgate_options = malloc( sizeof(int)*(nqubits) );

    // H on all qubits
    all_H = qdd_create_single_qubit_gates_same(nqubits, GATEID_H);
    qdd_protect(&all_H);

    // H on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_H;
    gate_list[n] = GATEID_I;
    first_n_H = qdd_create_single_qubit_gates(nqubits, gate_list);
    qdd_protect(&first_n_H);

    // oracle
    for (BDDVAR k = 0; k < n; k++) cgate_options[k] = flag[k];
    cgate_options[n] = 2; // ancilla qubit is target qubit of this multi-cgate
    oracle = qdd_create_multi_cgate(nqubits, cgate_options, GATEID_X);
    qdd_protect(&oracle);

    // X on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_X;
    gate_list[n] = GATEID_I;
    first_n_X = qdd_create_single_qubit_gates(nqubits, gate_list);
    qdd_protect(&first_n_X);

    // CZ gate on all qubits except ancilla
    for (BDDVAR k = 0; k < n-1; k++) cgate_options[k] = 1; // controls"
    cgate_options[n-1] = 2; // target last qubit before ancilla
    cgate_options[n] = -1; // do nothing to ancilla
    first_n_CZ = qdd_create_multi_cgate(nqubits, cgate_options, GATEID_Z);
    qdd_protect(&first_n_CZ);
    
    // Grover iteration = oracle + mean inversion
    grov_it = oracle;
    qdd_protect(&grov_it);
    grov_it = qdd_matmat_mult(first_n_H,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_CZ, grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = qdd_matmat_mult(first_n_H,  grov_it, nqubits);

    // Compute grov_it^t (by squaring if t is a power of 2)
    grov_its = grov_it;
    qdd_protect(&grov_its);
    if (EXP_BY_SQUARING && ((t & (t - 1)) == 0)) {
        for (int i = 0; i < log2l(t); i++) {
            grov_its = qdd_matmat_mult(grov_its, grov_its, nqubits);
        }
    }
    else {
        for (int i = 1; i < t; i++) {
            grov_its = qdd_matmat_mult(grov_its, grov_it, nqubits);
        }
    }
    R = R / t;
    *matrix = grov_its; // return for nodecount in bench

    qdd_unprotect(&first_n_H);
    qdd_unprotect(&first_n_X);
    qdd_unprotect(&first_n_CZ);
    qdd_unprotect(&oracle);
    qdd_unprotect(&grov_it);

    // Now, actually apply the circuit:
    // 1. Start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(nqubits) );
    for (BDDVAR k = 0; k < n; k++) init[k] = 0;
    init[n] = 1; 
    state = qdd_create_basis_state(nqubits, init);
    qdd_protect(&state);

    // 2. H on all qubits
    state = qdd_matvec_mult(all_H, state, nqubits);
    qdd_unprotect(&all_H);

    // 3. Grover iterations
    if (VERBOSE) {
        printf("\nGrover progress: %lf%%", 0.);
        fflush(stdout);
    }
    for (uint32_t i = 1; i <= R; i++) {
        if (i % 1000 == 0) {
            if (VERBOSE) {
                printf("\rGrover progress %lf%%", (double) i / (double) R * 100);
                fflush(stdout);
            }
            
            // gc amp table
            qdd_gc_amp_table();

            // gc node table
            sylvan_gc();
        }
        state = qdd_matvec_mult(grov_its, state, nqubits);
    }
    if (VERBOSE) {
        printf("\rGrover progress %lf%%\n", 100.);
    }

    qdd_unprotect(&grov_its);
    qdd_unprotect(&state);

    free(init);
    free(gate_list);
    free(cgate_options);

    return state;
}
