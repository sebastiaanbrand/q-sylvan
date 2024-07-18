#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include <sylvan_edge_weights_complex.h>
#include "test_assert.h"

bool VERBOSE = true;

// don't use for heap / malloced arrays
#define len(x) (sizeof(x) / sizeof(x[0]))

int test_x_gate()
{
    BDDVAR nqubits;
    QMDD v0, v1, v2, v3, v4, v5, mI, mX, mSqrtX, mSqrtXdag, mX1, mX0, mXXI, mTemp;
    bool x[] = {0};
    bool x3[] = {0, 0, 0};

    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qmdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qmdd_create_basis_state(nqubits, x);
    x[0] = 0; v2 = qmdd_create_basis_state(nqubits, x);
    mI = qmdd_create_all_identity_matrix(nqubits);
    mX = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_X);
    mSqrtX = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_sqrtX);
    mSqrtXdag = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_sqrtXdag);

    // matrix-vector mult
    v0 = aadd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v1);
    v0 = aadd_matvec_mult(mX, v0, nqubits); test_assert(v0 == v2);

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mI, mX, nqubits); test_assert(mTemp == mX);
    mTemp = aadd_matmat_mult(mX, mI, nqubits); test_assert(mTemp == mX);
    mTemp = aadd_matmat_mult(mX, mX, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mSqrtX,    mSqrtX,    nqubits); test_assert(mTemp == mX);
    mTemp = aadd_matmat_mult(mSqrtXdag, mSqrtXdag, nqubits); test_assert(mTemp == mX);
    mTemp = aadd_matmat_mult(mSqrtX,    mSqrtXdag, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mSqrtXdag, mSqrtX,    nqubits); test_assert(mTemp == mI);


    // 3 qubit test
    nqubits = 3;
    x3[2] = 0; x3[1] = 0; x3[0] = 0; v3 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; v4 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; v5 = qmdd_create_basis_state(nqubits, x3);
    mI  = qmdd_create_all_identity_matrix(nqubits);
    mX1 = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_X); // X gate on q1
    mX0 = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_X); // X gate on q0
    uint32_t gateids[] = {GATEID_X, GATEID_X, GATEID_I}; // X on q0 and q1, I on q2
    mXXI = qmdd_create_single_qubit_gates(nqubits, gateids);
    test_assert(aadd_countnodes(v3) == 4);
    test_assert(aadd_countnodes(v4) == 4);
    test_assert(aadd_countnodes(v5) == 4);
    
    // matrix-vector mult
    v3 = aadd_matvec_mult(mX1, v3, nqubits); test_assert(v3 == v4);
    v3 = aadd_matvec_mult(mX0, v3, nqubits); test_assert(v3 == v5);
    test_assert(aadd_countnodes(v3) == 4);
    test_assert(aadd_countnodes(v4) == 4);
    test_assert(aadd_countnodes(v5) == 4);

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mI, mX0, nqubits); test_assert(mTemp == mX0);
    mTemp = aadd_matmat_mult(mI, mX1, nqubits); test_assert(mTemp == mX1);
    mTemp = aadd_matmat_mult(mX0, mI, nqubits); test_assert(mTemp == mX0);
    mTemp = aadd_matmat_mult(mX1, mI, nqubits); test_assert(mTemp == mX1);
    mTemp = aadd_matmat_mult(mX0,mX0, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mX1,mX1, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mX0,mX1, nqubits); test_assert(mTemp == mXXI);
    mTemp = aadd_matmat_mult(mX1,mX0, nqubits); test_assert(mTemp == mXXI);

    // calculate (X0 X1)|00> by multiplying X0 and X1 first
    mTemp = aadd_matmat_mult(mX0, mX1, nqubits);
    v3 = qmdd_create_all_zero_state(nqubits);
    v3 = aadd_matvec_mult(mTemp, v3, nqubits);
    test_assert(v3 == v5);
    test_assert(aadd_countnodes(v5) == 4);


    if(VERBOSE) printf("matrix qmdd x gates:         ok\n");
    return 0;
}

int test_h_gate()
{
    BDDVAR nqubits;
    QMDD v0, v1, v2, v3, v4, v5, v6, mI, mH, mX, mZ, mII, mHI, mIH, mHH, mTemp;
    bool x[] = {0};
    bool x2[] = {0,0};
    AMP a;


    // Single qubit test
    nqubits = 1;
    x[0] = 0; v0 = qmdd_create_basis_state(nqubits, x);
    x[0] = 1; v1 = qmdd_create_basis_state(nqubits, x);
    mI = qmdd_create_all_identity_matrix(nqubits);
    mH = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mX = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_X);
    mZ = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_Z);

    // matrix-vector mult
    v0 = aadd_matvec_mult(mH, v0, nqubits);
    v1 = aadd_matvec_mult(mH, v1, nqubits);

    x[0] = 1; a = aadd_getvalue(v0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 0; a = aadd_getvalue(v0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 0; a = aadd_getvalue(v1, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 1; a = aadd_getvalue(v1, x); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mI, mH, nqubits); test_assert(mTemp == mH);
    mTemp = aadd_matmat_mult(mH, mI, nqubits); test_assert(mTemp == mH);
    mTemp = aadd_matmat_mult(mH, mH, nqubits); test_assert(mTemp == mI);

    // HXH = Z
    mTemp = mH;
    mTemp = aadd_matmat_mult(mTemp, mX, nqubits);
    mTemp = aadd_matmat_mult(mTemp, mH, nqubits);
    test_assert(mTemp == mZ);

    // HZH = X
    mTemp = mH;
    mTemp = aadd_matmat_mult(mTemp, mZ, nqubits);
    mTemp = aadd_matmat_mult(mTemp, mH, nqubits);
    test_assert(mTemp == mX);

    // Two qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v2 = qmdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 1; v3 = qmdd_create_basis_state(nqubits, x2); // |01>
    x2[1] = 0; x2[0] = 0; v4 = qmdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 0; v5 = qmdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 1; x2[0] = 0; v6 = qmdd_create_basis_state(nqubits, x2); // |10>
    mII = qmdd_create_all_identity_matrix(nqubits);
    mHI = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mIH = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mHH = qmdd_create_single_qubit_gates_same(nqubits, GATEID_H);

    // matrix-vector mult
    v2 = aadd_matvec_mult(mHI, v2, nqubits); // v2 = |0+>
    v3 = aadd_matvec_mult(mHI, v3, nqubits); // v3 = |0->
    v4 = aadd_matvec_mult(mIH, v4, nqubits); // v4 = |+0>
    v5 = aadd_matvec_mult(mHI, v5, nqubits);
    v5 = aadd_matvec_mult(mIH, v5, nqubits); // v5 = |++>
    v6 = aadd_matvec_mult(mHH, v6, nqubits); // v6 = |-+>

    // v2 = |0+>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v2, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v2, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(v2) == 2);

    // v3 = |0->
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v3, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v3, x2); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v3, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v3, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(v3) == 3);

    // v4 = |+0>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v4, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v4, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(v4) == 2);

    // v5 = |++>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v5, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(v5) == 1);

    // v6 = |-+>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v6, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v6, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v6, x2); test_assert(a == complex_lookup(-0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v6, x2); test_assert(a == complex_lookup(-0.5, 0));
    test_assert(aadd_countnodes(v6) == 2);

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mII, mHI, nqubits); test_assert(mTemp == mHI);
    mTemp = aadd_matmat_mult(mII, mIH, nqubits); test_assert(mTemp == mIH);
    mTemp = aadd_matmat_mult(mHI, mII, nqubits); test_assert(mTemp == mHI);
    mTemp = aadd_matmat_mult(mIH, mII, nqubits); test_assert(mTemp == mIH);
    mTemp = aadd_matmat_mult(mHI, mHI, nqubits); test_assert(mTemp == mII);
    mTemp = aadd_matmat_mult(mIH, mIH, nqubits); test_assert(mTemp == mII);
    mTemp = aadd_matmat_mult(mHI, mIH, nqubits); test_assert(mTemp == mHH);
    mTemp = aadd_matmat_mult(mIH, mHI, nqubits); test_assert(mTemp == mHH);
    mTemp = aadd_matmat_mult(mHH, mHH, nqubits); test_assert(mTemp == mII);

    // calculate (H0 H1)|00> by multiplying H0 and H1 first
    mTemp = aadd_matmat_mult(mHI, mIH, nqubits);
    v6 = qmdd_create_all_zero_state(nqubits);
    v6 = aadd_matvec_mult(mTemp, v6, nqubits);
    test_assert(v6 == v5);


    if(VERBOSE) printf("matrix qmdd h gates:         ok\n");
    return 0;
}

int test_phase_gates()
{
    BDDVAR nqubits;
    QMDD v0, vZ, vS, vSS, vT, vTT, vTTTT, vTTdag, vTdagT;
    QMDD mI, mH0, mH1, mZ0, mZ1, mS0, mSdag0, mT0, mT1, mTdag0, mTdag1, mTemp;
    bool x2[] = {0, 0};
    AMP a;


    // simple 2 qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qmdd_create_basis_state(nqubits, x2);
    mI     = qmdd_create_all_identity_matrix(nqubits);
    mH0    = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1    = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mZ0    = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_Z);
    mZ1    = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_Z);
    mS0    = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_S);
    mSdag0 = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_Sdag);
    mT0    = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_T);
    mT1    = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_T);
    mTdag0 = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_Tdag);
    mTdag1 = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_Tdag);
    v0 = aadd_matvec_mult(mH0, v0, nqubits);
    v0 = aadd_matvec_mult(mH1, v0, nqubits); // start with v0 = |++>

    // matrix-vector mult
    vZ     = aadd_matvec_mult(mZ0, v0, nqubits);
    vS     = aadd_matvec_mult(mS0, v0, nqubits);
    vSS    = aadd_matvec_mult(mS0, vS, nqubits);
    vT     = aadd_matvec_mult(mT0, v0, nqubits);
    vTT    = aadd_matvec_mult(mT0, vT, nqubits);
    vTTTT  = aadd_matvec_mult(mT0, vTT, nqubits);
    vTTTT  = aadd_matvec_mult(mT0, vTTTT, nqubits);
    vTTdag = aadd_matvec_mult(mT0, v0, nqubits);
    vTTdag = aadd_matvec_mult(mTdag0, vTTdag, nqubits);
    vTdagT = aadd_matvec_mult(mTdag0, v0, nqubits);
    vTdagT = aadd_matvec_mult(mT0, vTdagT, nqubits);
    test_assert(vZ == vSS);
    test_assert(vS == vTT);
    test_assert(vZ == vTTTT);
    test_assert(v0 == vTTdag);
    test_assert(v0 == vTdagT);

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mI,     mZ0,    nqubits); test_assert(mTemp == mZ0);
    mTemp = aadd_matmat_mult(mZ0,    mI,     nqubits); test_assert(mTemp == mZ0);
    mTemp = aadd_matmat_mult(mI,     mZ1,    nqubits); test_assert(mTemp == mZ1);
    mTemp = aadd_matmat_mult(mZ1,    mI,     nqubits); test_assert(mTemp == mZ1);
    mTemp = aadd_matmat_mult(mI,     mS0,    nqubits); test_assert(mTemp == mS0);
    mTemp = aadd_matmat_mult(mS0,    mI,     nqubits); test_assert(mTemp == mS0);
    mTemp = aadd_matmat_mult(mI,     mT0,    nqubits); test_assert(mTemp == mT0);
    mTemp = aadd_matmat_mult(mT0,    mI,     nqubits); test_assert(mTemp == mT0);
    mTemp = aadd_matmat_mult(mS0,    mS0,    nqubits); test_assert(mTemp == mZ0);
    mTemp = aadd_matmat_mult(mS0,    mSdag0, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mSdag0, mS0,    nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mSdag0, mSdag0, nqubits); test_assert(mTemp == mZ0);
    mTemp = aadd_matmat_mult(mT0,    mT0,    nqubits); test_assert(mTemp == mS0);
    mTemp = aadd_matmat_mult(mT1,    mTdag1, nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mTdag1, mT1,    nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mTdag0, mTdag0, nqubits); test_assert(mTemp == mSdag0);

    // T^7 == Tdag
    mTemp = qmdd_create_all_identity_matrix(nqubits);
    for (int k = 0; k < 7; k++) 
        mTemp = aadd_matmat_mult(mTemp, mT0, nqubits);
    test_assert(mTemp == mTdag0);


    // more matrix-vector mult    
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(v0) == 1);

    v0 = aadd_matvec_mult(mZ0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 2);

    v0 = aadd_matvec_mult(mZ0, v0, nqubits);
    v0 = aadd_matvec_mult(mZ1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 2);

    v0 = aadd_matvec_mult(mZ1, v0, nqubits);
    v0 = aadd_matvec_mult(mS0, v0, nqubits);
    v0 = aadd_matvec_mult(mS0, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 2);

    v0 = aadd_matvec_mult(mZ0, v0, nqubits);
    v0 = aadd_matvec_mult(mT1, v0, nqubits);
    v0 = aadd_matvec_mult(mT1, v0, nqubits);
    v0 = aadd_matvec_mult(mT1, v0, nqubits);
    v0 = aadd_matvec_mult(mT1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 2);

    v0 = aadd_matvec_mult(mZ1,    v0, nqubits);
    v0 = aadd_matvec_mult(mTdag1, v0, nqubits);
    v0 = aadd_matvec_mult(mTdag1, v0, nqubits);
    v0 = aadd_matvec_mult(mTdag1, v0, nqubits);
    v0 = aadd_matvec_mult(mTdag1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 2);
    

    if(VERBOSE) printf("matrix qmdd phase gates:     ok\n");
    return 0;
}

int test_cx_gate()
{
    BDDVAR nqubits;
    QMDD v0, mI, mH0, mH1, mCNOT, mCZ, mTemp;
    bool x2[] = {0,0};
    AMP a;


    // Test Bell state
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qmdd_create_basis_state(nqubits, x2);
    mI    = qmdd_create_all_identity_matrix(nqubits);
    mH0   = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1   = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mCNOT = qmdd_create_cgate(nqubits, 0, 1, GATEID_X);
    mCZ   = qmdd_create_cgate(nqubits, 0, 1, GATEID_Z);

    // matrix-matrix mult
    mTemp = aadd_matmat_mult(mI,    mH0,   nqubits); test_assert(mTemp == mH0);
    mTemp = aadd_matmat_mult(mH0,   mH0,   nqubits); test_assert(mTemp == mI);
    mTemp = aadd_matmat_mult(mI,    mCNOT, nqubits); test_assert(mTemp == mCNOT);
    mTemp = aadd_matmat_mult(mCNOT, mI,    nqubits); test_assert(mTemp == mCNOT);
    mTemp = aadd_matmat_mult(mCNOT, mCNOT, nqubits); test_assert(mTemp == mI);

    // H1 CNOT(0,1) H1 = CZ(0,1)
    mTemp = mH1;
    mTemp = aadd_matmat_mult(mTemp, mCNOT, nqubits);
    mTemp = aadd_matmat_mult(mTemp, mH1, nqubits);
    test_assert(mTemp == mCZ);

    // H1 CZ(0,1) H1 = CNOT(0,1)
    mTemp = mH1;
    mTemp = aadd_matmat_mult(mTemp, mCZ, nqubits);
    mTemp = aadd_matmat_mult(mTemp, mH1, nqubits);
    test_assert(mTemp == mCNOT);

    // matrix-vector mult
    v0 = aadd_matvec_mult(mH0, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(v0) == 2);

    v0 = aadd_matvec_mult(mCNOT, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    test_assert(aadd_countnodes(v0) == 4);

    // same as above but multiplies H CNOT first before applying to the state
    // note that we apply the H first, so it is on the right: CNOT H0 |00>
    mTemp = aadd_matmat_mult(mCNOT, mH0, nqubits);
    v0 = qmdd_create_all_zero_state(nqubits);
    v0 = aadd_matvec_mult(mTemp, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    test_assert(aadd_countnodes(v0) == 4);

    // If we did H0 CNOT |00> instead then the CNOT would not have an effect, 
    // and we'd just be left with H on qubit 0. This is also an example where
    // computing the circuit matrix first does a bunch of extra work if we
    // apply it on a state which is unaffected by (some of) the circuit.
    mTemp = aadd_matmat_mult(mH0, mCNOT, nqubits);
    v0 = qmdd_create_all_zero_state(nqubits);
    v0 = aadd_matvec_mult(mTemp, v0, nqubits);
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(v0) == 2);


    // TODO: test with slightly more qubits

    if(VERBOSE) printf("matrix qmdd cnot gates:      ok\n");
    return 0;
}

int test_cz_gate()
{
    BDDVAR nqubits;
    QMDD v0, mH0, mH1, mCZ;
    bool x2[] = {0, 0};
    AMP a;


    // 2 qubit graph state
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; v0 = qmdd_create_basis_state(nqubits, x2);
    mH0 = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1 = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mCZ = qmdd_create_cgate(nqubits, 0, 1, GATEID_Z);

    v0 = aadd_matvec_mult(mH0, v0, nqubits);
    v0 = aadd_matvec_mult(mH1, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(v0) == 1);

    v0 = aadd_matvec_mult(mCZ, v0, nqubits);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(v0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(v0) == 3);

    // TODO: test with slightly more qubits

    if(VERBOSE) printf("matrix qmdd cz gates:        ok\n");
    return 0;
}

int test_ccz_gate()
{
    BDDVAR nqubits;
    QMDD v3, vTemp, mCCZ, mH0, mH1, mH2;
    bool x3[] = {0,0,0};
    AMP a, aRef;


    // 3 qubit test
    nqubits = 3;
    v3 = qmdd_create_basis_state(nqubits, x3);
    mH0  = qmdd_create_single_qubit_gate(nqubits, 0, GATEID_H);
    mH1  = qmdd_create_single_qubit_gate(nqubits, 1, GATEID_H);
    mH2  = qmdd_create_single_qubit_gate(nqubits, 2, GATEID_H);

    v3 = aadd_matvec_mult(mH0, v3, nqubits);
    v3 = aadd_matvec_mult(mH1, v3, nqubits);
    v3 = aadd_matvec_mult(mH2, v3, nqubits);
    aRef = aadd_getvalue(v3, x3);
    test_assert(qmdd_is_unitvector(v3, nqubits));

    x3[2]=1; x3[1]=1; x3[0]=1;
    mCCZ = qmdd_create_all_control_phase(nqubits, x3);
    vTemp = aadd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == wgt_neg(aRef));

    x3[2]=0; x3[1]=1; x3[0]=1;
    mCCZ = qmdd_create_all_control_phase(nqubits, x3);
    vTemp = aadd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == wgt_neg(aRef));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);

    x3[2]=0; x3[1]=1; x3[0]=0;
    mCCZ = qmdd_create_all_control_phase(nqubits, x3);
    vTemp = aadd_matvec_mult(mCCZ, v3, nqubits);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == wgt_neg(aRef));
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(vTemp, x3); test_assert(a == aRef);


    if(VERBOSE) printf("matrix qmdd all-control z:   ok\n");
    return 0;
}

int test_multi_cgate()
{

    BDDVAR nqubits;
    QMDD qTest, qRef, qInit, matrix;


    uint32_t test_gates[] = {GATEID_X, GATEID_H, GATEID_Z, GATEID_sqrtX};

    // just single qubit gates
    nqubits = 3;
    qInit   = qmdd_create_all_zero_state(nqubits);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // t = 1
        int c_options[] = {-1,2,-1};
        qRef   = qmdd_gate(qInit, test_gates[i], 1);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest  = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // t = 0
        int c_options[] = {2,-1,-1};
        qRef   = qmdd_gate(qInit, test_gates[i], 0);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest  = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // single control gates
    nqubits = 5;
    qInit   = qmdd_create_all_zero_state(nqubits);
    qInit   = qmdd_gate(qInit, GATEID_H, 1);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c = 1, t = 2
        int c_options[] = {-1, 1, 2, -1, -1};
        qRef = qmdd_cgate(qInit, test_gates[i], 1, 2);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c = 1, t = 3
        int c_options[] = {-1, 1, -1, 2, -1};
        qRef = qmdd_cgate(qInit, test_gates[i], 1, 3);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // double control gates
    nqubits = 6;
    qInit   = qmdd_create_all_zero_state(nqubits);
    qInit   = qmdd_gate(qInit, GATEID_H, 0);
    qInit   = qmdd_gate(qInit, GATEID_H, 2);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // c1 = 0, c2 = 2, t = 5
        int c_options[] = {1, -1, 1, -1, -1, 2};
        qRef = qmdd_cgate2(qInit, test_gates[i], 0, 2, 5);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }
    for (uint32_t i = 0; i < len(test_gates); i++) { 
        // c1 = 0, c2 = 2, t = 5 (but now control c2 on |0> instead of |1>)
        int c_options[] = {1, -1, 0, -1, -1, 2};
        qRef = qmdd_gate(qInit, GATEID_X, 2);
        qRef = qmdd_cgate2(qRef, test_gates[i], 0, 2, 5);
        qRef = qmdd_gate(qRef, GATEID_X, 2);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // controls below target
    nqubits = 3;
    qInit   = qmdd_create_all_zero_state(nqubits);
    qInit   = qmdd_gate(qInit, GATEID_H, 0);
    qInit   = qmdd_gate(qInit, GATEID_X, 1);
    qInit   = qmdd_gate(qInit, GATEID_H, 2);
    for (uint32_t i = 0; i < len(test_gates); i++) {
        // t = 0, c1 = 1, c2 = 2 
        int c_options[] = {2, 1, 1};
        qRef = qmdd_circuit_swap(qInit, 0, 2);
        qRef = qmdd_cgate2(qRef, test_gates[i], 0, 1, 2);
        qRef = qmdd_circuit_swap(qRef, 0, 2);
        matrix = qmdd_create_multi_cgate(nqubits, c_options, test_gates[i]);
        qTest = aadd_matvec_mult(matrix, qInit, nqubits);
        test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
        test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
        test_assert(qTest == qRef);
    }

    // TODO: more tests

    if(VERBOSE) printf("matrix qmdd multi-cgate:     ok so far (WIP)\n");
    return 0;
}

int test_tensor_product()
{ 

    QMDD mTest, mRef, mX, mY, mZ, mH;

    uint32_t gateid_XY[]  = {GATEID_X, GATEID_Y};
    uint32_t gateid_XYZ[] = {GATEID_X, GATEID_Y, GATEID_Z};

    mX = qmdd_create_single_qubit_gate(1, 0, GATEID_X);
    mY = qmdd_create_single_qubit_gate(1, 0, GATEID_Y);
    mZ = qmdd_create_single_qubit_gate(1, 0, GATEID_Z);
    mH = qmdd_create_single_qubit_gate(1, 0, GATEID_H);

    // H (x) H
    mRef  = qmdd_create_single_qubit_gates_same(2, GATEID_H);
    mTest = aadd_mat_tensor_prod(mH, mH, 1);
    test_assert(aadd_is_ordered(mTest, 4)); // 2*n because of primed and unprimed
    test_assert(mTest == mRef);

    // X (x) Y
    mRef  = qmdd_create_single_qubit_gates(2, gateid_XY);
    mTest = aadd_mat_tensor_prod(mX, mY, 1);
    test_assert(aadd_is_ordered(mTest, 4));
    test_assert(mTest == mRef);

    // X  (x)  (Y (x) Z)
    mRef  = qmdd_create_single_qubit_gates(3, gateid_XYZ);
    mTest = aadd_mat_tensor_prod(mY, mZ, 1);
    mTest = aadd_mat_tensor_prod(mX, mTest, 1);
    test_assert(aadd_is_ordered(mTest, 6));
    test_assert(mTest == mRef);

    // (X (x) Y)  (x)  Z
    mRef  = qmdd_create_single_qubit_gates(3, gateid_XYZ);
    mTest = aadd_mat_tensor_prod(mX, mY, 1);
    mTest = aadd_mat_tensor_prod(mTest, mZ, 2);
    test_assert(aadd_is_ordered(mTest, 6));
    test_assert(mTest == mRef);

    // TODO: tests with bigger matrices

    if(VERBOSE) printf("matrix qmdd tensor product:  ok\n");
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

int test_with(int wgt_backend, int norm_strat, int wgt_indx_bits) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    double tolerance = -1; // default
    if (norm_strat == NORM_L2) tolerance = 1e-13;
    qsylvan_init_simulator(1LL<<wgt_indx_bits, 1LL<<wgt_indx_bits, tolerance,
                           wgt_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests
    aadd_set_auto_gc_wgt_table(false); // no auto gc of ctable yet for mult operations

    printf("wgt backend = %d, norm strat = %d, wgt indx bits = %d:\n", 
            wgt_backend, norm_strat, wgt_indx_bits);
    int res = runtests();

    sylvan_quit();
    lace_stop();

    return res;
}

int main()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
            if (test_with(backend, norm_strat, 11)) return 1;
            if (backend == COMP_HASHMAP) {
                // test with edge wgt index > 23 bits
                if (test_with(backend, norm_strat, 24)) return 1;
            }
        }
    }
    return 0;
}
