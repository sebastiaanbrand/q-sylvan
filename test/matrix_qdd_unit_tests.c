#include <stdio.h>
#include <time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_complex.h"

bool VERBOSE = true;

// don't use for heap / malloced arrays
#define len(x) (sizeof(x) / sizeof(x[0]))

int test_x_gate()
{
    BDDVAR nqubits;
    QDD v0, v1, v2, v3, v4, v5, mI, mX, mX1, mX0, mXXI, mTemp;
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
    mTemp = qdd_matmat_mult(mX, mI, nqubits); test_assert(mTemp == mX);
    mTemp = qdd_matmat_mult(mX, mX, nqubits); test_assert(mTemp == mI);


    // 3 qubit test
    nqubits = 3;
    x3[2] = 0; x3[1] = 0; x3[0] = 0; v3 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; v4 = qdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; v5 = qdd_create_basis_state(nqubits, x3);
    mI  = qdd_create_all_identity_matrix(nqubits);
    mX1 = qdd_create_single_qubit_gate(nqubits, 1, GATEID_X); // X gate on q1
    mX0 = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X); // X gate on q0
    uint32_t gateids[] = {GATEID_X, GATEID_X, GATEID_I}; // X on q0 and q1, I on q2
    mXXI = qdd_create_single_qubit_gates(nqubits, gateids);
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
    mTemp = qdd_matmat_mult(mX0, mI, nqubits); test_assert(mTemp == mX0);
    mTemp = qdd_matmat_mult(mX1, mI, nqubits); test_assert(mTemp == mX1);
    mTemp = qdd_matmat_mult(mX0,mX0, nqubits); test_assert(mTemp == mI);
    mTemp = qdd_matmat_mult(mX1,mX1, nqubits); test_assert(mTemp == mI);
    mTemp = qdd_matmat_mult(mX0,mX1, nqubits); test_assert(mTemp == mXXI);
    mTemp = qdd_matmat_mult(mX1,mX0, nqubits); test_assert(mTemp == mXXI);

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
    QDD v0, v1, v2, v3, v4, v5, v6, mI, mH, mX, mZ, mII, mHI, mIH, mHH, mTemp;
    bool x[] = {0};
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qdd_create_basis_state(nqubits, x);
    mI = qdd_create_all_identity_matrix(nqubits);
    mH = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mX = qdd_create_single_qubit_gate(nqubits, 0, GATEID_X);
    mZ = qdd_create_single_qubit_gate(nqubits, 0, GATEID_Z);

    // matrix-vector mult
    v0 = qdd_matvec_mult(mH, v0, nqubits);
    v1 = qdd_matvec_mult(mH, v1, nqubits);

    x[0] = 1; a = qdd_get_amplitude(v0, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 0; a = qdd_get_amplitude(v0, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 0; a = qdd_get_amplitude(v1, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 1; a = qdd_get_amplitude(v1, x); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mI, mH, nqubits); test_assert(mTemp == mH);
    mTemp = qdd_matmat_mult(mH, mI, nqubits); test_assert(mTemp == mH);
    mTemp = qdd_matmat_mult(mH, mH, nqubits); test_assert(mTemp == mI);

    // HXH = Z
    mTemp = mH;
    mTemp = qdd_matmat_mult(mTemp, mX, nqubits);
    mTemp = qdd_matmat_mult(mTemp, mH, nqubits);
    test_assert(mTemp == mZ);

    // HZH = X
    mTemp = mH;
    mTemp = qdd_matmat_mult(mTemp, mZ, nqubits);
    mTemp = qdd_matmat_mult(mTemp, mH, nqubits);
    test_assert(mTemp == mX);

    // Two qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v2 = qdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 1; v3 = qdd_create_basis_state(nqubits, x2); // |01>
    x2[1] = 0; x2[0] = 0; v4 = qdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 0; v5 = qdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 1; x2[0] = 0; v6 = qdd_create_basis_state(nqubits, x2); // |10>
    mII = qdd_create_all_identity_matrix(nqubits);
    mHI = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mIH = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mHH = qdd_create_single_qubit_gates_same(nqubits, GATEID_H);

    // matrix-vector mult
    v2 = qdd_matvec_mult(mHI, v2, nqubits); // v2 = |0+>
    v3 = qdd_matvec_mult(mHI, v3, nqubits); // v3 = |0->
    v4 = qdd_matvec_mult(mIH, v4, nqubits); // v4 = |+0>
    v5 = qdd_matvec_mult(mHI, v5, nqubits);
    v5 = qdd_matvec_mult(mIH, v5, nqubits); // v5 = |++>
    v6 = qdd_matvec_mult(mHH, v6, nqubits); // v6 = |-+>

    // v2 = |0+>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v2, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v2, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v2, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v2, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v2) == 2);

    // v3 = |0->
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v3, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v3, x2); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v3, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v3, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v3) == 3);

    // v4 = |+0>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v4, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v4, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v4, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v4, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v4) == 2);

    // v5 = |++>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(v5) == 1);

    // v6 = |-+>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v6, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v6, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v6, x2); test_assert(a == comp_lookup(comp_make(-0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v6, x2); test_assert(a == comp_lookup(comp_make(-0.5, 0)));
    test_assert(qdd_countnodes(v6) == 2);

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mII, mHI, nqubits); test_assert(mTemp == mHI);
    mTemp = qdd_matmat_mult(mII, mIH, nqubits); test_assert(mTemp == mIH);
    mTemp = qdd_matmat_mult(mHI, mII, nqubits); test_assert(mTemp == mHI);
    mTemp = qdd_matmat_mult(mIH, mII, nqubits); test_assert(mTemp == mIH);
    mTemp = qdd_matmat_mult(mHI, mHI, nqubits); test_assert(mTemp == mII);
    mTemp = qdd_matmat_mult(mIH, mIH, nqubits); test_assert(mTemp == mII);
    mTemp = qdd_matmat_mult(mHI, mIH, nqubits); test_assert(mTemp == mHH);
    mTemp = qdd_matmat_mult(mIH, mHI, nqubits); test_assert(mTemp == mHH);
    mTemp = qdd_matmat_mult(mHH, mHH, nqubits); test_assert(mTemp == mII);

    // calculate (H0 H1)|00> by multiplying H0 and H1 first
    mTemp = qdd_matmat_mult(mHI, mIH, nqubits);
    v6 = qdd_create_all_zero_state(nqubits);
    v6 = qdd_matvec_mult(mTemp, v6, nqubits);
    test_assert(v6 == v5);


    if(VERBOSE) printf("matrix qdd h gates:         ok\n");
    return 0;
}

int test_phase_gates()
{
    BDDVAR nqubits;
    QDD v0, vZ, vS, vSS, vT, vTT, vTTTT, vTTdag, vTdagT;
    QDD mI, mH0, mH1, mZ0, mZ1, mS0, mT0, mT1, mTdag0, mTdag1, mTemp;
    bool x2[] = {0, 0};
    AMP a;

    LACE_ME;

    // simple 2 qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qdd_create_basis_state(nqubits, x2);
    mI     = qdd_create_all_identity_matrix(nqubits);
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

    // matrix-vector mult
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

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mI,     mZ0,    nqubits); test_assert(mTemp == mZ0);
    mTemp = qdd_matmat_mult(mZ0,    mI,     nqubits); test_assert(mTemp == mZ0);
    mTemp = qdd_matmat_mult(mI,     mZ1,    nqubits); test_assert(mTemp == mZ1);
    mTemp = qdd_matmat_mult(mZ1,    mI,     nqubits); test_assert(mTemp == mZ1);
    mTemp = qdd_matmat_mult(mI,     mS0,    nqubits); test_assert(mTemp == mS0);
    mTemp = qdd_matmat_mult(mS0,    mI,     nqubits); test_assert(mTemp == mS0);
    mTemp = qdd_matmat_mult(mI,     mT0,    nqubits); test_assert(mTemp == mT0);
    mTemp = qdd_matmat_mult(mT0,    mI,     nqubits); test_assert(mTemp == mT0);
    mTemp = qdd_matmat_mult(mS0,    mS0,    nqubits); test_assert(mTemp == mZ0);
    mTemp = qdd_matmat_mult(mT0,    mT0,    nqubits); test_assert(mTemp == mS0);
    mTemp = qdd_matmat_mult(mT1,    mTdag1, nqubits); test_assert(mTemp == mI);
    mTemp = qdd_matmat_mult(mTdag1, mT1,    nqubits); test_assert(mTemp == mI);

    // T^7 == Tdag
    mTemp = qdd_create_all_identity_matrix(nqubits);
    for (int k = 0; k < 7; k++) 
        mTemp = qdd_matmat_mult(mTemp, mT0, nqubits);
    test_assert(mTemp == mTdag0);


    // more matrix-vector mult    
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(v0) == 1);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);
    v0 = qdd_matvec_mult(mZ1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ1, v0, nqubits);
    v0 = qdd_matvec_mult(mS0, v0, nqubits);
    v0 = qdd_matvec_mult(mS0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ0, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);
    v0 = qdd_matvec_mult(mT1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mZ1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);
    v0 = qdd_matvec_mult(mTdag1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 2);
    

    if(VERBOSE) printf("matrix qdd phase gates:     ok\n");
    return 0;
}

int test_cx_gate()
{
    BDDVAR nqubits;
    QDD v0, mI, mH0, mH1, mCNOT, mCZ, mTemp;
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Test Bell state
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qdd_create_basis_state(nqubits, x2);
    mI    = qdd_create_all_identity_matrix(nqubits);
    mH0   = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1   = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mCNOT = qdd_create_controlled_gate(nqubits, 0, 1, GATEID_X);
    mCZ   = qdd_create_controlled_gate(nqubits, 0, 1, GATEID_Z);

    // matrix-matrix mult
    mTemp = qdd_matmat_mult(mI,    mH0,   nqubits); test_assert(mTemp == mH0);
    mTemp = qdd_matmat_mult(mH0,   mH0,   nqubits); test_assert(mTemp == mI);
    mTemp = qdd_matmat_mult(mI,    mCNOT, nqubits); test_assert(mTemp == mCNOT);
    mTemp = qdd_matmat_mult(mCNOT, mI,    nqubits); test_assert(mTemp == mCNOT);
    mTemp = qdd_matmat_mult(mCNOT, mCNOT, nqubits); test_assert(mTemp == mI);

    // H1 CNOT(0,1) H1 = CZ(0,1)
    mTemp = mH1;
    mTemp = qdd_matmat_mult(mTemp, mCNOT, nqubits);
    mTemp = qdd_matmat_mult(mTemp, mH1, nqubits);
    test_assert(mTemp == mCZ);

    // H1 CZ(0,1) H1 = CNOT(0,1)
    mTemp = mH1;
    mTemp = qdd_matmat_mult(mTemp, mCZ, nqubits);
    mTemp = qdd_matmat_mult(mTemp, mH1, nqubits);
    test_assert(mTemp == mCNOT);

    // matrix-vector mult
    v0 = qdd_matvec_mult(mH0, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v0) == 2);

    v0 = qdd_matvec_mult(mCNOT, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    test_assert(qdd_countnodes(v0) == 4);

    // same as above but multiplies H CNOT first before applying to the state
    // note that we apply the H first, so it is on the right: CNOT H0 |00>
    mTemp = qdd_matmat_mult(mCNOT, mH0, nqubits);
    v0 = qdd_create_all_zero_state(nqubits);
    v0 = qdd_matvec_mult(mTemp, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    test_assert(qdd_countnodes(v0) == 4);

    // If we did H0 CNOT |00> instead then the CNOT would not have an effect, 
    // and we'd just be left with H on qubit 0. This is also an example where
    // computing the circuit matrix first does a bunch of extra work if we
    // apply it on a state which is unaffected by (some of) the circuit.
    mTemp = qdd_matmat_mult(mH0, mCNOT, nqubits);
    v0 = qdd_create_all_zero_state(nqubits);
    v0 = qdd_matvec_mult(mTemp, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(v0) == 2);


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

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(v0) == 1);

    v0 = qdd_matvec_mult(mCZ, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(v0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(v0) == 3);

    // TODO: test with slightly more qubits

    if(VERBOSE) printf("matrix qdd cz gates:        ok\n");
    return 0;
}

int test_ccz_gate()
{
    BDDVAR nqubits;
    QDD v3, vTemp, mCCZ, mH0, mH1, mH2;
    bool x3[] = {0,0,0};
    AMP a, aRef;

    LACE_ME;

    // 3 qubit test
    nqubits = 3;
    v3 = qdd_create_basis_state(nqubits, x3);
    mH0  = qdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1  = qdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mH2  = qdd_create_single_qubit_gate(nqubits, 2, GATEID_H);

    v3 = qdd_matvec_mult(mH0, v3, nqubits);
    v3 = qdd_matvec_mult(mH1, v3, nqubits);
    v3 = qdd_matvec_mult(mH2, v3, nqubits);
    aRef = qdd_get_amplitude(v3, x3);
    test_assert(qdd_is_unitvector(v3, nqubits));

    x3[2]=1; x3[1]=1; x3[0]=1;
    mCCZ = qdd_create_all_control_phase(nqubits, x3);
    vTemp = qdd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == amp_neg(aRef));

    x3[2]=0; x3[1]=1; x3[0]=1;
    mCCZ = qdd_create_all_control_phase(nqubits, x3);
    vTemp = qdd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == amp_neg(aRef));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);

    x3[2]=0; x3[1]=1; x3[0]=0;
    mCCZ = qdd_create_all_control_phase(nqubits, x3);
    vTemp = qdd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == amp_neg(aRef));
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(vTemp, x3); test_assert(a == aRef);


    if(VERBOSE) printf("matrix qdd all-control z:   ok\n");
    return 0;
}

int test_multi_cgate()
{

    BDDVAR nqubits;
    QDD qTest, qRef, qInit, matrix;

    LACE_ME;

    uint32_t test_gates[] = {GATEID_X, GATEID_H, GATEID_Z, GATEID_sqrtX};

    // just single qubit gates
    nqubits = 3;
    qInit   = qdd_create_all_zero_state(nqubits);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // t = 1
        int c_options[] = {-1,2,-1};
        qRef   = qdd_gate(qInit, test_gates[i], 1);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest  = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // t = 0
        int c_options[] = {2,-1,-1};
        qRef   = qdd_gate(qInit, test_gates[i], 0);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest  = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // single control gates
    nqubits = 5;
    qInit   = qdd_create_all_zero_state(nqubits);
    qInit   = qdd_gate(qInit, GATEID_H, 1);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c = 1, t = 2
        int c_options[] = {-1, 1, 2, -1, -1};
        qRef = qdd_cgate(qInit, test_gates[i], 1, 2);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c = 1, t = 3
        int c_options[] = {-1, 1, -1, 2, -1};
        qRef = qdd_cgate(qInit, test_gates[i], 1, 3);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // double control gates
    nqubits = 6;
    qInit   = qdd_create_all_zero_state(nqubits);
    qInit   = qdd_gate(qInit, GATEID_H, 0);
    qInit   = qdd_gate(qInit, GATEID_H, 2);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c1 = 0, c2 = 2, t = 5
        int c_options[] = {1, -1, 1, -1, -1, 2};
        qRef = qdd_cgate2(qInit, test_gates[i], 0, 2, 5);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) { 
        // c1 = 0, c2 = 2, t = 5 (but now control c2 on |0> instead of |1>)
        int c_options[] = {1, -1, 0, -1, -1, 2};
        qRef = qdd_gate(qInit, GATEID_X, 2);
        qRef = qdd_cgate2(qRef, test_gates[i], 0, 2, 5);
        qRef = qdd_gate(qRef, GATEID_X, 2);
        matrix = qdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = qdd_matvec_mult(matrix, qInit, nqubits);
        test_assert(qdd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(qdd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }


    // TODO: more tests


    if(VERBOSE) printf("matrix qdd multi-cgate:     ok so far (WIP)\n");
    return 0;
}

int test_tensor_product()
{
    LACE_ME; 

    QDD mTest, mRef, mX, mY, mZ, mH;

    uint32_t gateid_XY[]  = {GATEID_X, GATEID_Y};
    uint32_t gateid_XYZ[] = {GATEID_X, GATEID_Y, GATEID_Z};

    mX = qdd_create_single_qubit_gate(1, 0, GATEID_X);
    mY = qdd_create_single_qubit_gate(1, 0, GATEID_Y);
    mZ = qdd_create_single_qubit_gate(1, 0, GATEID_Z);
    mH = qdd_create_single_qubit_gate(1, 0, GATEID_H);

    // H (x) H
    mRef  = qdd_create_single_qubit_gates_same(2, GATEID_H);
    mTest = qdd_mat_tensor_prod(mH, mH, 1);
    test_assert(qdd_is_ordered(mTest, 4)); // 2*n because of primed and unprimed
    test_assert(mTest == mRef);

    // X (x) Y
    mRef  = qdd_create_single_qubit_gates(2, gateid_XY);
    mTest = qdd_mat_tensor_prod(mX, mY, 1);
    test_assert(qdd_is_ordered(mTest, 4));
    test_assert(mTest == mRef);

    // X  (x)  (Y (x) Z)
    mRef  = qdd_create_single_qubit_gates(3, gateid_XYZ);
    mTest = qdd_mat_tensor_prod(mY, mZ, 1);
    mTest = qdd_mat_tensor_prod(mX, mTest, 1);
    test_assert(qdd_is_ordered(mTest, 6));
    test_assert(mTest == mRef);

    // (X (x) Y)  (x)  Z
    mRef  = qdd_create_single_qubit_gates(3, gateid_XYZ);
    mTest = qdd_mat_tensor_prod(mX, mY, 1);
    mTest = qdd_mat_tensor_prod(mTest, mZ, 2);
    test_assert(qdd_is_ordered(mTest, 6));
    test_assert(mTest == mRef);

    // TODO: tests with bigger matrices

    if(VERBOSE) printf("matrix qdd tensor product:  ok\n");
    return 0;
}

int runtests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    if (test_x_gate()) return 1;
    if (test_h_gate()) return 1;
    if (test_phase_gates()) return 1;
    if (test_cx_gate()) return 1;
    if (test_cz_gate()) return 1;
    if (test_ccz_gate()) return 1;
    if (test_multi_cgate()) return 1;
    if (test_tensor_product()) return 1;

    return 0;
}

int test_with_cmap()
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<11, -1, COMP_HASHMAP);
    qdd_set_testing_mode(true); // turn on internal sanity tests
    qdd_set_auto_gc_amp_table(false); // no auto gc of ctable yet for mult operations

    printf("using cmap:\n");
    int res = runtests();

    sylvan_quit();
    lace_exit();

    return res;
}

int test_with_rmap()
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<11, -1, REAL_HASHMAP);
    qdd_set_testing_mode(true); // turn on internal sanity tests
    qdd_set_auto_gc_amp_table(false); // no auto gc of ctable yet for mult operations

    printf("using rmap:\n");
    int res = runtests();

    sylvan_quit();
    lace_exit();

    return res;
}

int test_with_tree_map()
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<11, -1, REAL_TREE);
    qdd_set_testing_mode(true); // turn on internal sanity tests
    qdd_set_auto_gc_amp_table(false); // no auto gc of ctable yet for mult operations

    printf("using tree map:\n");
    int res = runtests();

    sylvan_quit();
    lace_exit();

    return res;
}

int main()
{
    if (test_with_cmap()) return 1;
    if (test_with_rmap()) return 1;
    if (test_with_tree_map()) return 1;
    return 0;
}
