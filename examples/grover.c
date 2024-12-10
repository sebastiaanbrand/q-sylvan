#include "grover.h"

static bool VERBOSE = false;
static bool EXP_BY_SQUARING = true;

void
qmdd_grover_set_verbose(bool v)
{
    VERBOSE = v;
}

uint64_t
qmdd_grover_approx_number_of_gates(BDDVAR nbits)
{
    uint64_t iterations = floor( 3.14159265359/4.0 * sqrt( pow(2,nbits) ) );
    uint64_t init_gates = (nbits + 1);
    uint64_t gts_per_it = (4 + 4*nbits);
    return (init_gates + iterations*gts_per_it);
}

bool *
qmdd_grover_random_flag(BDDVAR nbits)
{
    bool *flag = malloc( (sizeof(bool)*nbits) );
    for (BDDVAR k = 0; k < nbits; k++) flag[k] = (rand() % 2);
    return flag;
}

bool *
qmdd_grover_ones_flag(BDDVAR nbits)
{
    bool *flag = malloc( (sizeof(bool)*nbits) );
    for (BDDVAR k = 0; k < nbits; k++) flag[k] = 1;
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


TASK_IMPL_3(QMDD, qmdd_grover_iteration, QMDD, qmdd, BDDVAR, n, bool*, oracle)
{
    // "oracle" call  (apply -1 flag to desired amplitude)
    qmdd = qmdd_gate(qmdd, GATEID_H, n);
    qmdd = qmdd_all_control_phase(qmdd, n+1, oracle);
    qmdd = qmdd_gate(qmdd, GATEID_H, n);

    // H on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qmdd = qmdd_gate(qmdd, GATEID_H, k);

    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qmdd = qmdd_gate(qmdd, GATEID_X, k);

    // Controlled Z over all qubits (except ancilla)
    qmdd = qmdd_cgate_range(qmdd, GATEID_Z, 0, n-2, n-1);

    // X on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qmdd = qmdd_gate(qmdd, GATEID_X, k);

    // H on all qubits (except ancilla)
    for (BDDVAR k = 0; k < n; k++) qmdd = qmdd_gate(qmdd, GATEID_H, k);
    return qmdd;
}

QMDD
qmdd_grover(BDDVAR n, bool *flag)
{
    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(n+1) );
    for (BDDVAR k = 0; k < n; k++) init[k] = 0;
    init[n] = 1; 
    QMDD qmdd = qmdd_create_basis_state(n+1, init);
    free(init);

    // H on all qubits
    for (BDDVAR k = 0; k < n+1; k++) qmdd = qmdd_gate(qmdd, GATEID_H, k);

    // Grover iterations
    bool *oracle = append_one(n, flag);
    for (uint32_t i = 1; i <= R; i++) {
        qmdd = qmdd_grover_iteration(qmdd, n, oracle);
    }
    free(oracle);

    return qmdd;
}

QMDD
qmdd_grover_matrix(BDDVAR n, bool *flag)
{
    QMDD m;
    return qmdd_grover_matrix_multi_its(n, flag, 1, &m);
}

QMDD qmdd_grover_matrix_multi_its(BDDVAR n, bool *flag, int t, QMDD *matrix)
{
    BDDVAR nqubits = n + 1;

    // Number of Grover iterations
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // Create matrix QMDDs
    QMDD all_H, first_n_H, first_n_X, oracle, first_n_CZ, grov_it, grov_its, state;
    uint32_t *gate_list = malloc( sizeof(uint32_t)*(nqubits) );
    int *cgate_options = malloc( sizeof(int)*(nqubits) );

    // H on all qubits
    all_H = qmdd_create_single_qubit_gates_same(nqubits, GATEID_H);
    evbdd_protect(&all_H);

    // H on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_H;
    gate_list[n] = GATEID_I;
    first_n_H = qmdd_create_single_qubit_gates(nqubits, gate_list);
    evbdd_protect(&first_n_H);

    // oracle
    for (BDDVAR k = 0; k < n; k++) cgate_options[k] = flag[k];
    cgate_options[n] = 2; // ancilla qubit is target qubit of this multi-cgate
    oracle = qmdd_create_multi_cgate(nqubits, cgate_options, GATEID_X);
    evbdd_protect(&oracle);

    // X on first n qubits (all qubits except ancilla)
    for (BDDVAR k = 0; k < n; k++) gate_list[k] = GATEID_X;
    gate_list[n] = GATEID_I;
    first_n_X = qmdd_create_single_qubit_gates(nqubits, gate_list);
    evbdd_protect(&first_n_X);

    // CZ gate on all qubits except ancilla
    for (BDDVAR k = 0; k < n-1; k++) cgate_options[k] = 1; // controls"
    cgate_options[n-1] = 2; // target last qubit before ancilla
    cgate_options[n] = -1; // do nothing to ancilla
    first_n_CZ = qmdd_create_multi_cgate(nqubits, cgate_options, GATEID_Z);
    evbdd_protect(&first_n_CZ);
    
    // Grover iteration = oracle + mean inversion
    grov_it = oracle;
    evbdd_protect(&grov_it);
    grov_it = evbdd_matmat_mult(first_n_H,  grov_it, nqubits);
    grov_it = evbdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = evbdd_matmat_mult(first_n_CZ, grov_it, nqubits);
    grov_it = evbdd_matmat_mult(first_n_X,  grov_it, nqubits);
    grov_it = evbdd_matmat_mult(first_n_H,  grov_it, nqubits);

    // Compute grov_it^t (by squaring if t is a power of 2)
    grov_its = grov_it;
    evbdd_protect(&grov_its);
    if (EXP_BY_SQUARING && ((t & (t - 1)) == 0)) {
        for (int i = 0; i < log2l(t); i++) {
            grov_its = evbdd_matmat_mult(grov_its, grov_its, nqubits);
        }
    }
    else {
        for (int i = 1; i < t; i++) {
            grov_its = evbdd_matmat_mult(grov_its, grov_it, nqubits);
        }
    }
    R = R / t;
    *matrix = grov_its; // return for nodecount in bench

    evbdd_unprotect(&first_n_H);
    evbdd_unprotect(&first_n_X);
    evbdd_unprotect(&first_n_CZ);
    evbdd_unprotect(&oracle);
    evbdd_unprotect(&grov_it);

    // Now, actually apply the circuit:
    // 1. Start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(nqubits) );
    for (BDDVAR k = 0; k < n; k++) init[k] = 0;
    init[n] = 1; 
    state = qmdd_create_basis_state(nqubits, init);
    evbdd_protect(&state);

    // 2. H on all qubits
    state = evbdd_matvec_mult(all_H, state, nqubits);
    evbdd_unprotect(&all_H);

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
            
            // gc edge wegith table
            evbdd_gc_wgt_table();

            // gc node table
            sylvan_gc();
        }
        state = evbdd_matvec_mult(grov_its, state, nqubits);
    }
    if (VERBOSE) {
        printf("\rGrover progress %lf%%\n", 100.);
    }

    evbdd_unprotect(&grov_its);
    evbdd_unprotect(&state);

    free(init);
    free(gate_list);
    free(cgate_options);

    return state;
}
