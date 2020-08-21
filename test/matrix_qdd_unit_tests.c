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
    QDD v0, v1, v2, v3, v4, v5, mI, mX, mX1, mX0, mTemp;
    bool x[] = {0};
    bool x3[] = {0, 0, 0};

    LACE_ME;

    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qdd_create_basis_state(nqubits, x);
    x[0] = 0; v2 = qdd_create_basis_state(nqubits, x);
    mI = qdd_create_all_identity_matrix(nqubits);
    mX = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X);

    // matrix-vector mult
    v0 = qdd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v1);
    v0 = qdd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v2);

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mI, mX, nqubits); test_assert(mTemp == mX);
    mTemp = qdd_matmat_mult(mX, mX, nqubits); test_assert(mTemp == mI);

    // 3 qubit test
    nqubits = 3;
    x3[2] = 0; x3[1] = 0; x3[0] = 0; v3 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; v4 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; v5 = qdd_create_basis_state(nqubits, x3);
    mI  = qdd_create_all_identity_matrix(nqubits);
    mX1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_X); // X gate on q1
    mX0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X); // X gate on q0
    test_assert(qdd_countnodes(v3) == 4);
    test_assert(qdd_countnodes(v4) == 4);
    test_assert(qdd_countnodes(v5) == 4);
    
    // matrix-vector mult
    v3 = qdd_matvec_mult(mX1, v3, nqubits); test_assert(v3 == v4);
    v3 = qdd_matvec_mult(mX0, v3, nqubits); test_assert(v3 == v5);
    test_assert(qdd_countnodes(v3) == 4);
    test_assert(qdd_countnodes(v4) == 4);
    test_assert(qdd_countnodes(v5) == 4);

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mI, mX0, nqubits); test_assert(mTemp == mX0);
    mTemp = qdd_matmat_mult(mI, mX1, nqubits); test_assert(mTemp == mX1);
    mTemp = qdd_matmat_mult(mX0,mX0, nqubits); test_assert(mTemp == mI);
    mTemp = qdd_matmat_mult(mX1,mX1, nqubits); test_assert(mTemp == mI);

    // calculate (X0 X1)|00> by multiplying X0 and X1 first
    mTemp = qdd_matmat_mult(mX0,mX1, nqubits);
    v3 = qdd_create_all_zero_state(nqubits);
    v3 = qdd_matvec_mult(mTemp, v3, nqubits);
    test_assert(v3 == v5);
    test_assert(qdd_countnodes(v5) == 4);


    if(VERBOSE) printf("matrix qdd x gates:         ok\n");
    return 0;
}

int test_h_gate()
{
    BDDVAR nqubits;
    QDD v0, v1, v2, v3, v4, v5, mH, mH0, mH1;
    bool x[] = {0};
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qdd_create_basis_state(nqubits, x);
    mH = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);

    v0 = qdd_matvec_mult(mH, v0, nqubits);
    v1 = qdd_matvec_mult(mH, v1, nqubits);

    x[0] = 1; a = qdd_get_amplitude(v0, x); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x[0] = 0; a = qdd_get_amplitude(v0, x); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x[0] = 0; a = qdd_get_amplitude(v1, x); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x[0] = 1; a = qdd_get_amplitude(v1, x); test_assert(a == Clookup(Cmake(Qmake(0,-1,2),0)));


    // Two qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v2 = qdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 1; v3 = qdd_create_basis_state(nqubits, x2); // |01>
    x2[1] = 0; x2[0] = 0; v4 = qdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 0; v5 = qdd_create_basis_state(nqubits, x2); // |00>
    mH0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);

    v2 = qdd_matvec_mult(mH0, v2, nqubits); // v2 = |0+>
    v3 = qdd_matvec_mult(mH0, v3, nqubits); // v3 = |0->
    v4 = qdd_matvec_mult(mH1, v4, nqubits); // v4 = |+0>
    v5 = qdd_matvec_mult(mH0, v5, nqubits);
    v5 = qdd_matvec_mult(mH1, v5, nqubits); // v5 = |++>

    // v2 = |0+>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v2, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v2, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v2, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v2, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v2) == 2);

    // v3 = |0->
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v3, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v3, x2); test_assert(a == Clookup(Cmake(Qmake(0,-1,2),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v3, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v3, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v3) == 3);

    // v4 = |+0>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v4, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v4, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v4, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v4, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v4) == 2);

    // v5 = |++>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v5, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v5, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v5, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v5, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    test_assert(qdd_countnodes(v5) == 1);


    if(VERBOSE) printf("matrix qdd h gates:         ok\n");
    return 0;
}

int test_phase_gates()
{
    BDDVAR nqubits;
    QDD v0, vZ, vS, vSS, vT, vTT, vTTTT, vTTdag, vTdagT;
    QDD mH0, mH1, mZ0, mZ1, mS0, mT0, mT1, mTdag0, mTdag1;
    bool x2[] = {0, 0};
    AMP a;

    LACE_ME;

    // simple 2 qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qdd_create_basis_state(nqubits, x2);
    mH0    = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1    = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mZ0    = qdd_create_single_qubit_gate(nqubits, 0, GATEID_Z);
    mZ1    = qdd_create_single_qubit_gate(nqubits, 1, GATEID_Z);
    mS0    = qdd_create_single_qubit_gate(nqubits, 0, GATEID_S);
    mT0    = qdd_create_single_qubit_gate(nqubits, 0, GATEID_T);
    mT1    = qdd_create_single_qubit_gate(nqubits, 1, GATEID_T);
    mTdag0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_Tdag);
    mTdag1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_Tdag);
    v0 = qdd_matvec_mult(mH0, v0, nqubits);
    v0 = qdd_matvec_mult(mH1, v0, nqubits); // start with v0 = |++>

    vZ     = qdd_matvec_mult(mZ0, v0, nqubits);
    vS     = qdd_matvec_mult(mS0, v0, nqubits);
    vSS    = qdd_matvec_mult(mS0, vS, nqubits);
    vT     = qdd_matvec_mult(mT0, v0, nqubits);
    vTT    = qdd_matvec_mult(mT0, vT, nqubits);
    vTTTT  = qdd_matvec_mult(mT0, vTT, nqubits);
    vTTTT  = qdd_matvec_mult(mT0, vTTTT, nqubits);
    vTTdag = qdd_matvec_mult(mT0, v0, nqubits);
    vTTdag = qdd_matvec_mult(mTdag0, vTTdag, nqubits);
    vTdagT = qdd_matvec_mult(mTdag0, v0, nqubits);
    vTdagT = qdd_matvec_mult(mT0, vTdagT, nqubits);

    test_assert(vZ == vSS);
    test_assert(vS == vTT);
    test_assert(vZ == vTTTT);
    test_assert(v0 == vTTdag);
    test_assert(v0 == vTdagT);

    
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    test_assert(qdd_countnodes(v0) == 1);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);
    v0 = qdd_matvec_mult(mZ1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ1, v0, nqubits);
    v0 = qdd_matvec_mult(mS0, v0, nqubits);
    v0 = qdd_matvec_mult(mS0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits)

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);
    

    if(VERBOSE) printf("matrix qdd phase gates:     ok\n");
    return 0;
}

int test_cx_gate()
{
    BDDVAR nqubits;
    QDD v0, mH, mCNOT;
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Test Bell state
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qdd_create_basis_state(nqubits, x2);
    mH    = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mCNOT = qdd_create_controlled_gate(nqubits, 0, 1, GATEID_X);

    v0 = qdd_matvec_mult(mH, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mCNOT, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(Qmake(0,1,2),0)));
    test_assert(qdd_countnodes(v0) == 4);

    // TODO: test with slightly more qubits

    if(VERBOSE) printf("matrix qdd cnot gates:      ok\n");
    return 0;
}

int test_cz_gate()
{
    BDDVAR nqubits;
    QDD v0, mH0, mH1, mCZ;
    bool x2[] = {0, 0};
    AMP a;

    LACE_ME;

    // 2 qubit graph state
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qdd_create_basis_state(nqubits, x2);
    mH0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mCZ = qdd_create_controlled_gate(nqubits, 0, 1, GATEID_Z);

    v0 = qdd_matvec_mult(mH0, v0, nqubits);
    v0 = qdd_matvec_mult(mH1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    test_assert(qdd_countnodes(v0) == 1);

    v0 = qdd_matvec_mult(mCZ, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == Clookup(Cmake(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 3);

    // TODO: test with slightly more qubits

    if(VERBOSE) printf("matrix qdd cz gates:        ok\n");
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
