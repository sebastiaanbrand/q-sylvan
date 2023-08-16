#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"
#include <sylvan_edge_weights_complex.h>

bool VERBOSE = true;

int test_x_gate()
{
    QMDD q0, q1, q2, q3, q4, q5;
    bool x[] = {0};
    bool x3[] = {0, 0, 0};
    BDDVAR nqubits;

    // Single qubit test
    x[0] = 0; q0 = qmdd_create_basis_state(1, x);
    x[0] = 1; q1 = qmdd_create_basis_state(1, x);
    x[0] = 0; q2 = qmdd_create_basis_state(1, x);

    q0 = qmdd_gate(q0, GATEID_X, 0); test_assert(q0 == q1);
    q0 = qmdd_gate(q0, GATEID_X, 0); test_assert(q0 == q2);
    q0 = qmdd_gate(q0, GATEID_X, 0); test_assert(q0 == q1);
    q0 = qmdd_gate(q0, GATEID_X, 0); test_assert(q0 == q2);

    // 3 qubit test
    nqubits = 3;
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; q4 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; q5 = qmdd_create_basis_state(nqubits, x3);
    test_assert(aadd_countnodes(q3) == 4);
    test_assert(aadd_countnodes(q4) == 4);
    test_assert(aadd_countnodes(q5) == 4);
    test_assert(aadd_is_ordered(q3, nqubits));
    test_assert(aadd_is_ordered(q4, nqubits));
    test_assert(aadd_is_ordered(q5, nqubits));
    
    q3 = qmdd_gate(q3, GATEID_X, 1); test_assert(q3 == q4);
    q3 = qmdd_gate(q3, GATEID_X, 0); test_assert(q3 == q5);
    test_assert(aadd_countnodes(q3) == 4);
    test_assert(aadd_countnodes(q4) == 4);
    test_assert(aadd_countnodes(q5) == 4);
    test_assert(aadd_is_ordered(q3, nqubits));
    test_assert(aadd_is_ordered(q4, nqubits));
    test_assert(aadd_is_ordered(q5, nqubits));
    
    // Same 3 qubit test with sqrt(X)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; q4 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; q5 = qmdd_create_basis_state(nqubits, x3);

    q3 = qmdd_gate(q3, GATEID_sqrtX, 1); test_assert(qmdd_is_unitvector(q3, nqubits));
    q3 = qmdd_gate(q3, GATEID_sqrtX, 1); test_assert(qmdd_is_unitvector(q3, nqubits));
    test_assert(q3 == q4);

    q3 = qmdd_gate(q3, GATEID_sqrtX, 0); test_assert(qmdd_is_unitvector(q3, nqubits));
    q3 = qmdd_gate(q3, GATEID_sqrtX, 0); test_assert(qmdd_is_unitvector(q3, nqubits));
    test_assert(q3 == q5);
    
    test_assert(aadd_countnodes(q3) == 4);
    test_assert(aadd_countnodes(q4) == 4);
    test_assert(aadd_countnodes(q5) == 4);
    test_assert(aadd_is_ordered(q3, nqubits));
    test_assert(aadd_is_ordered(q4, nqubits));
    test_assert(aadd_is_ordered(q5, nqubits));

    // Same 3 qubit test with sqrt(X)^dag
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; q4 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; q5 = qmdd_create_basis_state(nqubits, x3);

    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 1); test_assert(qmdd_is_unitvector(q3, nqubits));
    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 1); test_assert(qmdd_is_unitvector(q3, nqubits));
    test_assert(q3 == q4);

    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 0); test_assert(qmdd_is_unitvector(q3, nqubits));
    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 0); test_assert(qmdd_is_unitvector(q3, nqubits));
    test_assert(q3 == q5);
    
    test_assert(aadd_countnodes(q3) == 4);
    test_assert(aadd_countnodes(q4) == 4);
    test_assert(aadd_countnodes(q5) == 4);
    test_assert(aadd_is_ordered(q3, nqubits));
    test_assert(aadd_is_ordered(q4, nqubits));
    test_assert(aadd_is_ordered(q5, nqubits));

    // Test sqrt(X) sqrt(X)^dag gives I
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qmdd_create_basis_state(nqubits, x3);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q4 = qmdd_create_basis_state(nqubits, x3);

    q3 = qmdd_gate(q3, GATEID_sqrtX, 1);
    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 1);
    test_assert(q3 == q4);

    q3 = qmdd_gate(q3, GATEID_sqrtXdag, 1);
    q3 = qmdd_gate(q3, GATEID_sqrtX, 1);
    test_assert(q3 == q4);


    if(VERBOSE) printf("qmdd x gates:              ok\n");
    return 0;
}

int test_h_gate()
{
    QMDD q0, q1, q2, q3, q4, q5;
    bool x[] = {0};
    bool x2[] = {0,0};
    AMP a;
    BDDVAR nqubits;

    // Single qubit test
    x[0] = 0; q0 = qmdd_create_basis_state(1, x);
    x[0] = 1; q1 = qmdd_create_basis_state(1, x);

    q0 = qmdd_gate(q0, GATEID_H, 0);
    q1 = qmdd_gate(q1, GATEID_H, 0);

    x[0] = 0; a = aadd_getvalue(q0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 1; a = aadd_getvalue(q0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 0; a = aadd_getvalue(q1, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 1; a = aadd_getvalue(q1, x); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));


    // Two qubit test
    nqubits = 2;
    x2[1] = 0; x2[0] = 0; q2 = qmdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 1; q3 = qmdd_create_basis_state(nqubits, x2); // |01>
    x2[1] = 0; x2[0] = 0; q4 = qmdd_create_basis_state(nqubits, x2); // |00>
    x2[1] = 0; x2[0] = 0; q5 = qmdd_create_basis_state(nqubits, x2); // |00>
    q2 = qmdd_gate(q2, GATEID_H, 0); // q2 = |0+>
    q3 = qmdd_gate(q3, GATEID_H, 0); // q3 = |0->
    q4 = qmdd_gate(q4, GATEID_H, 1); // q4 = |+0>
    q5 = qmdd_gate(q5, GATEID_H, 0);
    q5 = qmdd_gate(q5, GATEID_H, 1); // q5 = |++>
    test_assert(aadd_is_ordered(q2, nqubits));
    test_assert(aadd_is_ordered(q3, nqubits));
    test_assert(aadd_is_ordered(q4, nqubits));
    test_assert(aadd_is_ordered(q5, nqubits));

    // q2 = |0+>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q2, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q2, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(q2) == 2);

    // q3 = |0->
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q3, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q3, x2); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q3, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q3, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(q3) == 3);

    // q4 = |+0>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q4, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q4, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(q4) == 2);

    // q5 = |++>
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(q5) == 1);

    if(VERBOSE) printf("qmdd h gates:              ok\n");
    return 0;
}

int test_phase_gates()
{
    QMDD q0, qZ, qS, qSS, qSdag, qSdagS, qT, qTT, qTTTT, qTTdag, qTdagT;
    bool x2[] = {0, 0};
    AMP a;

    // simple 2 qubit test
    x2[1] = 0; x2[0] = 0; q0 = qmdd_create_basis_state(2, x2);
    q0 = qmdd_gate(q0, GATEID_H, 0);
    q0 = qmdd_gate(q0, GATEID_H, 1);

    qZ     = qmdd_gate(q0, GATEID_Z, 0);
    qS     = qmdd_gate(q0, GATEID_S, 0);
    qSdag  = qmdd_gate(q0, GATEID_Sdag, 0);
    qSdagS = qmdd_gate(qSdag, GATEID_S, 0);
    qSS    = qmdd_gate(qS, GATEID_S, 0);
    qT     = qmdd_gate(q0, GATEID_T, 0);
    qTT    = qmdd_gate(qT, GATEID_T, 0);
    qTTTT  = qmdd_gate(qTT, GATEID_T, 0);
    qTTTT  = qmdd_gate(qTTTT, GATEID_T, 0);
    qTTdag = qmdd_gate(q0, GATEID_T, 0);
    qTTdag = qmdd_gate(qTTdag, GATEID_Tdag, 0);
    qTdagT = qmdd_gate(q0, GATEID_Tdag, 0);
    qTdagT = qmdd_gate(qTdagT, GATEID_T, 0);

    test_assert(qZ == qSS);
    test_assert(q0 == qSdagS);
    test_assert(qS == qTT);
    test_assert(qZ == qTTTT);
    test_assert(q0 == qTTdag);
    test_assert(q0 == qTdagT);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(q0) == 1);

    q0 = qmdd_gate(q0, GATEID_Z, 0);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 0);
    q0 = qmdd_gate(q0, GATEID_Z, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_S, 0);
    q0 = qmdd_gate(q0, GATEID_S, 0);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 0);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_Sdag, 1);
    q0 = qmdd_gate(q0, GATEID_Sdag, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(q0) == 2);

    // check R_k gates
    test_assert(gates[GATEID_Rk(0)][3] == gates[GATEID_I][3]);
    test_assert(gates[GATEID_Rk(1)][3] == gates[GATEID_Z][3]);
    test_assert(gates[GATEID_Rk(2)][3] == gates[GATEID_S][3]);
    test_assert(gates[GATEID_Rk(3)][3] == gates[GATEID_T][3]);
    test_assert(gates[GATEID_Rk_dag(0)][3] == gates[GATEID_I][3]);
    test_assert(gates[GATEID_Rk_dag(1)][3] == gates[GATEID_Z][3]);
    test_assert(gates[GATEID_Rk_dag(2)][3] == gates[GATEID_Sdag][3]);
    test_assert(gates[GATEID_Rk_dag(3)][3] == gates[GATEID_Tdag][3]);

    if(VERBOSE) printf("qmdd phase gates:          ok\n");
    return 0;
}

int test_pauli_rotation_gates()
{
    QMDD qInit, qTest, qRef;
    BDDVAR nqubits, t;

    // Rz rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);
    qInit = qmdd_gate(qInit, GATEID_H, t);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(1.0), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // Z gate
    qRef  = qmdd_gate(qInit, GATEID_Z, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(0.5), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // S gate
    qRef  = qmdd_gate(qInit, GATEID_S, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(0.25), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // T gate
    qRef  = qmdd_gate(qInit, GATEID_T, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(0.125), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    // Rx rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(1.0), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // X gate
    qRef  = qmdd_gate(qInit, GATEID_X, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(0.5), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(X) gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtX, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(0.25), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    // Ry rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);
    qInit = qmdd_gate(qInit, GATEID_H, t);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(1.0), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // Y gate
    qRef  = qmdd_gate(qInit, GATEID_Y, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(0.5), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(Y) gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtY, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(0.25), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(Y)^dag gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtYdag, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(-0.25), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(aadd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(aadd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    if(VERBOSE) printf("qmdd Rx, Ry, Rz gates:     ok\n");
    return 0;
}

int test_cx_gate()
{
    QMDD qBell;
    bool x2[] = {0,0};
    AMP a;

    // Test Bell state
    x2[1] = 0; x2[0] = 0; qBell = qmdd_create_basis_state(2, x2);
    qBell = qmdd_gate(qBell, GATEID_H, 0);
    
    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(qBell, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(qBell, x2); test_assert(a == AADD_ZERO);
    test_assert(aadd_countnodes(qBell) == 2);

    qBell = qmdd_cgate(qBell, GATEID_X, 0, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(qBell, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(qBell, x2); test_assert(a == AADD_ZERO);
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    test_assert(aadd_countnodes(qBell) == 4);

    // TODO: more tests

    if(VERBOSE) printf("qmdd cnot gates:           ok\n");
    return 0;
}

int test_cz_gate()
{
    QMDD qGraph;
    bool x2[] = {0, 0};
    AMP a;

    // 2 qubit graph state
    x2[1] = 0; x2[0] = 0; qGraph =qmdd_create_basis_state(2, x2);
    qGraph = qmdd_gate(qGraph, GATEID_H, 0);
    qGraph = qmdd_gate(qGraph, GATEID_H, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(aadd_countnodes(qGraph) == 1);

    qGraph = qmdd_cgate(qGraph, GATEID_Z, 0, 1);

    x2[1] = 0; x2[0] = 0; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = aadd_getvalue(qGraph, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(aadd_countnodes(qGraph) == 3);

    if(VERBOSE) printf("qmdd CZ gates:             ok\n");
    return 0;
}

int test_ccz_gate()
{
    QMDD q3;
    bool x3[] = {0,0,0};
    AMP a, aRef;

    q3 = qmdd_create_basis_state(3, x3);
    q3 = qmdd_gate(q3, GATEID_H, 0);
    q3 = qmdd_gate(q3, GATEID_H, 1);
    q3 = qmdd_gate(q3, GATEID_H, 2);
    aRef = aadd_getvalue(q3, x3);

    complex_t mone = cmone();
    x3[2]=1; x3[1]=0; x3[0]=1; q3 = qmdd_all_control_phase(q3, 3, x3);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == aRef);    
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == wgt_mul(aRef,weight_lookup(&mone)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == aRef);
    test_assert(aadd_is_ordered(q3, 3));

    // TODO: few more tests

    if(VERBOSE) printf("qmdd all-control z gate:   ok\n");
    return 0;
}

int test_controlled_range_gate()
{
    QMDD q10, qres;
    bool x10[] = {0,0,0,0,0,0,0,0,0,0};
    AMP a, aRef, aRefMin;
    bool *x_bits;

    complex_t mone = cmone();
    BDDVAR nqubits = 10;
    q10 = qmdd_create_basis_state(10, x10);
    for (BDDVAR k = 0; k < nqubits; k++) q10 = qmdd_gate(q10, GATEID_H, k);
    aRef = aadd_getvalue(q10, x10);
    aRefMin = wgt_mul(aRef, weight_lookup(&mone));

    // assert |++..+> state
    test_assert(aadd_is_ordered(q10, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = aadd_getvalue(q10, x_bits); 
        test_assert(a == aRef);
        free(x_bits); // int_to_bitarray mallocs
    }

    // CZ gate on c=0,1,2,3 t=7
    qres = qmdd_cgate_range(q10, GATEID_Z, 0, 3, 7);
    test_assert(aadd_is_ordered(qres, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = aadd_getvalue(qres, x_bits);
        // for all amps |1111***1**> we should have a phase flip
        if (x_bits[0] && x_bits[1] && x_bits[2] && x_bits[3] && x_bits[7]) {
            test_assert(a == aRefMin);
        }
        // all others should remain the same
        else {
            test_assert(a == aRef);
        }
        free(x_bits); // int_to_bitarray mallocs
    }

    // CZ gate on c=2,3,4,5 t=8
    qres = qmdd_cgate_range(q10, GATEID_Z, 2, 5, 8);
    test_assert(aadd_is_ordered(qres, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = aadd_getvalue(qres, x_bits);
        // for all amps |**1111**1*> we should have a phase flip
        if (x_bits[2] && x_bits[3] && x_bits[4] && x_bits[5] && x_bits[8]) {
            test_assert(a == aRefMin);
        }
        // all others should remain the same
        else {
            test_assert(a == aRef);
        }
        free(x_bits); // int_to_bitarray mallocs
    }

    // TODO: test other gates as well

    if(VERBOSE) printf("qmdd controlled range:     ok\n");
    return 0;
}

int run_qmdd_tests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    // gates
    if (test_x_gate()) return 1;
    if (test_h_gate()) return 1;
    if (test_phase_gates()) return 1;
    if (test_pauli_rotation_gates()) return 1;
    if (test_cx_gate()) return 1;
    if (test_cz_gate()) return 1;
    if (test_controlled_range_gate()) return 1;
    if (test_ccz_gate()) return 1;

    return 0;
}

int test_with(int amps_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<11, -1, amps_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("amps backend = %d, norm strategy = %d:\n", amps_backend, norm_strat);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
            if (test_with(backend, norm_strat)) return 1;
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
