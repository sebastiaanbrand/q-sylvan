#include <stdio.h>
#include <time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_int.h"

bool VERBOSE = true;

int test_matrix_creation()
{
    QDD matrix = qdd_create_single_qubit_gate(3, 1, GATEID_Z);
    FILE *fp;
    fp = fopen("matrix_test.dot", "w");
    qdd_fprintdot(fp, matrix, false);

    if(VERBOSE) printf("matrix creation:            TODO\n");
    return 0;
}

int test_x_gate()
{
    BDDVAR nqubits;
    QDD v0, v1, v2, v3, v4, v5, mX, mX1, mX0;
    bool x[] = {0};
    bool x3[] = {0, 0, 0};

    LACE_ME;

    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qdd_create_basis_state(nqubits, x);
    x[0] = 0; v2 = qdd_create_basis_state(nqubits, x);
    mX = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X);

    v0 = qdd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v1);
    v0 = qdd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v2);

    // 3 qubit test
    nqubits = 3;
    x3[2] = 0; x3[1] = 0; x3[0] = 0; v3 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; v4 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; v5 = qdd_create_basis_state(nqubits, x3);
    mX1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_X); // X gate on q1
    mX0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X); // X gate on q0
    test_assert(qdd_countnodes(v3) == 4);
    test_assert(qdd_countnodes(v4) == 4);
    test_assert(qdd_countnodes(v5) == 4);
    
    v3 = qdd_matvec_mult(mX1, v3, nqubits); test_assert(v3 == v4);
    v3 = qdd_matvec_mult(mX0, v3, nqubits); test_assert(v3 == v5);
    test_assert(qdd_countnodes(v3) == 4);
    test_assert(qdd_countnodes(v4) == 4);
    test_assert(qdd_countnodes(v5) == 4);

    if(VERBOSE) printf("matrix qdd x gates:         ok\n");
    return 0;
}

int test_h_gate()
{


    if(VERBOSE) printf("matrix qdd z gates:         TODO\n");
    return 0;
}

int test_phase_gates()
{


    if(VERBOSE) printf("matrix qdd phase gates:     TODO\n");
    return 0;
}

int test_cx_gate()
{


    if(VERBOSE) printf("matrix qdd cnot gates:      TODO\n");
    return 0;
}

int test_cz_gate()
{


    if(VERBOSE) printf("matrix qdd cz gates:        TODO\n");
    return 0;
}

int test_ccz_gate()
{


    if(VERBOSE) printf("matrix qdd all-control cz:  TODO\n");
    return 0;
}

int runtests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    if (test_matrix_creation()) return 1;
    if (test_x_gate()) return 1;
    if (test_h_gate()) return 1;
    if (test_phase_gates()) return 1;
    if (test_cx_gate()) return 1;
    if (test_cz_gate()) return 1;
    if (test_ccz_gate()) return 1;

    return 0;
}

int main()
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s)\n", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<19);

    int res = runtests();

    free_amplitude_table();
    sylvan_quit();
    lace_exit();

    return res;
}
