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
    test_assert(evbdd_countnodes(q3) == 4);
    test_assert(evbdd_countnodes(q4) == 4);
    test_assert(evbdd_countnodes(q5) == 4);
    test_assert(evbdd_is_ordered(q3, nqubits));
    test_assert(evbdd_is_ordered(q4, nqubits));
    test_assert(evbdd_is_ordered(q5, nqubits));
    
    q3 = qmdd_gate(q3, GATEID_X, 1); test_assert(q3 == q4);
    q3 = qmdd_gate(q3, GATEID_X, 0); test_assert(q3 == q5);
    test_assert(evbdd_countnodes(q3) == 4);
    test_assert(evbdd_countnodes(q4) == 4);
    test_assert(evbdd_countnodes(q5) == 4);
    test_assert(evbdd_is_ordered(q3, nqubits));
    test_assert(evbdd_is_ordered(q4, nqubits));
    test_assert(evbdd_is_ordered(q5, nqubits));
    
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
    
    test_assert(evbdd_countnodes(q3) == 4);
    test_assert(evbdd_countnodes(q4) == 4);
    test_assert(evbdd_countnodes(q5) == 4);
    test_assert(evbdd_is_ordered(q3, nqubits));
    test_assert(evbdd_is_ordered(q4, nqubits));
    test_assert(evbdd_is_ordered(q5, nqubits));

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
    
    test_assert(evbdd_countnodes(q3) == 4);
    test_assert(evbdd_countnodes(q4) == 4);
    test_assert(evbdd_countnodes(q5) == 4);
    test_assert(evbdd_is_ordered(q3, nqubits));
    test_assert(evbdd_is_ordered(q4, nqubits));
    test_assert(evbdd_is_ordered(q5, nqubits));

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

    x[0] = 0; a = evbdd_getvalue(q0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 1; a = evbdd_getvalue(q0, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 0; a = evbdd_getvalue(q1, x); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x[0] = 1; a = evbdd_getvalue(q1, x); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));


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
    test_assert(evbdd_is_ordered(q2, nqubits));
    test_assert(evbdd_is_ordered(q3, nqubits));
    test_assert(evbdd_is_ordered(q4, nqubits));
    test_assert(evbdd_is_ordered(q5, nqubits));

    // q2 = |0+>
    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q2, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q2, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q2, x2); test_assert(a == EVBDD_ZERO);
    test_assert(evbdd_countnodes(q2) == 2);

    // q3 = |0->
    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q3, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q3, x2); test_assert(a == complex_lookup(-1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q3, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q3, x2); test_assert(a == EVBDD_ZERO);
    test_assert(evbdd_countnodes(q3) == 3);

    // q4 = |+0>
    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q4, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q4, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q4, x2); test_assert(a == EVBDD_ZERO);
    test_assert(evbdd_countnodes(q4) == 2);

    // q5 = |++>
    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q5, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(evbdd_countnodes(q5) == 1);

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

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(evbdd_countnodes(q0) == 1);

    q0 = qmdd_gate(q0, GATEID_Z, 0);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 0);
    q0 = qmdd_gate(q0, GATEID_Z, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_S, 0);
    q0 = qmdd_gate(q0, GATEID_S, 0);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 0);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);
    q0 = qmdd_gate(q0, GATEID_T, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);
    q0 = qmdd_gate(q0, GATEID_Tdag, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

    q0 = qmdd_gate(q0, GATEID_Z, 1);
    q0 = qmdd_gate(q0, GATEID_Sdag, 1);
    q0 = qmdd_gate(q0, GATEID_Sdag, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(q0, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(q0) == 2);

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

    double pi = 2.0 * flt_acos(0.0);

    // Rz rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);
    qInit = qmdd_gate(qInit, GATEID_H, t);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(2.0*pi), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // Z gate
    qRef  = qmdd_gate(qInit, GATEID_Z, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(pi), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // S gate
    qRef  = qmdd_gate(qInit, GATEID_S, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(pi/2.0), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // T gate
    qRef  = qmdd_gate(qInit, GATEID_T, t);
    qTest = qmdd_gate(qInit, GATEID_Rz(pi/4.0), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    // Phase(theta) gate (similar to Rz(theta) but different global phase)
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);
    qInit = qmdd_gate(qInit, GATEID_H, t);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Phase(2.0*pi), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // Z gate
    qRef  = qmdd_gate(qInit, GATEID_Z, t);
    qTest = qmdd_gate(qInit, GATEID_Phase(pi), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // S gate
    qRef  = qmdd_gate(qInit, GATEID_S, t);
    qTest = qmdd_gate(qInit, GATEID_Phase(pi/2.0), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // T gate
    qRef  = qmdd_gate(qInit, GATEID_T, t);
    qTest = qmdd_gate(qInit, GATEID_Phase(pi/4.0), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    // Rx rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(2.0*pi), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // X gate
    qRef  = qmdd_gate(qInit, GATEID_X, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(pi), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(X) gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtX, t);
    qTest = qmdd_gate(qInit, GATEID_Rx(pi/2.0), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // U(theta,-pi/2,pi/2) == Rx(theta)
    qRef  = qmdd_gate(qInit, GATEID_Rx(pi/2.0), t);
    qTest = qmdd_gate(qInit, GATEID_U(pi/2.0, -pi/2.0, pi/2.0), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // U(theta,-pi/2,pi/2) == Rx(theta)
    qRef  = qmdd_gate(qInit, GATEID_Rx(1.42), t);
    qTest = qmdd_gate(qInit, GATEID_U(1.42, -pi/2.0, pi/2.0), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    // Ry rotations
    nqubits = 3, t = 1;
    qInit = qmdd_create_all_zero_state(nqubits);
    qInit = qmdd_gate(qInit, GATEID_H, t);

    // I gate
    qRef  = qmdd_gate(qInit, GATEID_I, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(2.0*pi), t);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // Y gate
    qRef  = qmdd_gate(qInit, GATEID_Y, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(pi), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(Y) gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtY, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(pi/2.0), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // sqrt(Y)^dag gate
    qRef  = qmdd_gate(qInit, GATEID_sqrtYdag, t);
    qTest = qmdd_gate(qInit, GATEID_Ry(-pi/2.0), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // U(theta,0,0) == Ry(theta)
    qRef  = qmdd_gate(qInit, GATEID_Ry(pi/2.0), t);
    qTest = qmdd_gate(qInit, GATEID_U(pi/2.0, 0, 0), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // U(theta,0,0) == Ry(theta)
    qRef  = qmdd_gate(qInit, GATEID_Ry(0.66), t);
    qTest = qmdd_gate(qInit, GATEID_U(0.66, 0, 0), t);
    qRef  = qmdd_remove_global_phase(qRef);
    qTest = qmdd_remove_global_phase(qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, false));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    if(VERBOSE) printf("qmdd Rx, Ry, Rz gates:     ok\n");
    return 0;
}

int test_cx_gate()
{
    QMDD qBell;
    BDDVAR nqubits = 2;
    bool x2[] = {0,0};
    AMP a;

    // Test Bell state
    x2[1] = 0; x2[0] = 0; qBell = qmdd_create_basis_state(nqubits, x2);
    qBell = qmdd_gate(qBell, GATEID_H, 0);
    
    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    test_assert(evbdd_countnodes(qBell) == 2);

    qBell = qmdd_cgate(qBell, GATEID_X, 0, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    test_assert(evbdd_countnodes(qBell) == 4);

    // Test Bell state with CX upside down
    x2[1] = 0; x2[0] = 0; qBell = qmdd_create_basis_state(nqubits, x2);
    qBell = qmdd_gate(qBell, GATEID_H, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);

    qBell = qmdd_cgate(qBell, GATEID_X, 1, 0, nqubits);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qBell, x2); test_assert(a == EVBDD_ZERO);
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qBell, x2); test_assert(a == complex_lookup(1.0/flt_sqrt(2.0),0));

    // TODO: more tests

    if(VERBOSE) printf("qmdd cnot gates:           ok\n");
    return 0;
}

int test_cz_gate()
{
    QMDD qGraph;
    BDDVAR nqubits = 2;
    bool x2[] = {0, 0};
    AMP a;

    // 2 qubit graph state
    x2[1] = 0; x2[0] = 0; qGraph =qmdd_create_basis_state(nqubits, x2);
    qGraph = qmdd_gate(qGraph, GATEID_H, 0);
    qGraph = qmdd_gate(qGraph, GATEID_H, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(evbdd_countnodes(qGraph) == 1);

    qGraph = qmdd_cgate(qGraph, GATEID_Z, 0, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(qGraph) == 3);

    // same 2 qubit graph state but with control and target arguments swapped
    x2[1] = 0; x2[0] = 0; qGraph =qmdd_create_basis_state(nqubits, x2);
    qGraph = qmdd_gate(qGraph, GATEID_H, 0);
    qGraph = qmdd_gate(qGraph, GATEID_H, 1);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    test_assert(evbdd_countnodes(qGraph) == 1);

    qGraph = qmdd_cgate(qGraph, GATEID_Z, 1, 0, nqubits);

    x2[1] = 0; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 0; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 0; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(0.5, 0));
    x2[1] = 1; x2[0] = 1; a = evbdd_getvalue(qGraph, x2); test_assert(a == complex_lookup(-0.5,0));
    test_assert(evbdd_countnodes(qGraph) == 3);

    if(VERBOSE) printf("qmdd CZ gates:             ok\n");
    return 0;
}

int test_ccz_gate()
{
    QMDD q3;
    BDDVAR nqubits = 3;
    bool x3[] = {0,0,0};
    AMP a, aRef;
    complex_t mone = cmone();

    // Test CCZ using qmdd_all_control_phase()
    q3 = qmdd_create_all_zero_state(nqubits);
    q3 = qmdd_gate(q3, GATEID_H, 0);
    q3 = qmdd_gate(q3, GATEID_H, 1);
    q3 = qmdd_gate(q3, GATEID_H, 2);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; aRef = evbdd_getvalue(q3, x3);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qmdd_all_control_phase(q3, nqubits, x3);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);    
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == wgt_mul(aRef,weight_lookup(&mone)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    test_assert(evbdd_is_ordered(q3, 3));

    // Test CCZ using qmdd_cgate2()
    q3 = qmdd_create_all_zero_state(nqubits);
    q3 = qmdd_gate(q3, GATEID_H, 0);
    q3 = qmdd_gate(q3, GATEID_H, 1);
    q3 = qmdd_gate(q3, GATEID_H, 2);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; aRef = evbdd_getvalue(q3, x3);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qmdd_cgate2(q3, GATEID_Z, 0, 1, 2);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);    
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == wgt_mul(aRef,weight_lookup(&mone)));
    test_assert(evbdd_is_ordered(q3, 3));

    // Test CCZ using qmdd_cgate2() with target between controls
    q3 = qmdd_create_all_zero_state(nqubits);
    q3 = qmdd_gate(q3, GATEID_H, 0);
    q3 = qmdd_gate(q3, GATEID_H, 1);
    q3 = qmdd_gate(q3, GATEID_H, 2);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; aRef = evbdd_getvalue(q3, x3);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qmdd_cgate2(q3, GATEID_Z, 0, 2, 1, nqubits);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);    
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == wgt_mul(aRef,weight_lookup(&mone)));
    test_assert(evbdd_is_ordered(q3, 3));

    // TODO: more tests?

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
    aRef = evbdd_getvalue(q10, x10);
    aRefMin = wgt_mul(aRef, weight_lookup(&mone));

    // assert |++..+> state
    test_assert(evbdd_is_ordered(q10, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = evbdd_getvalue(q10, x_bits); 
        test_assert(a == aRef);
        free(x_bits); // int_to_bitarray mallocs
    }

    // CZ gate on c=0,1,2,3 t=7
    qres = qmdd_cgate_range(q10, GATEID_Z, 0, 3, 7);
    test_assert(evbdd_is_ordered(qres, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = evbdd_getvalue(qres, x_bits);
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
    test_assert(evbdd_is_ordered(qres, nqubits));
    for (uint64_t x = 0; x < (1UL<<nqubits); x++) {
        x_bits = int_to_bitarray(x, nqubits, true);
        a = evbdd_getvalue(qres, x_bits);
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

int test_with(int wgt_backend, int norm_strat, int wgt_indx_bits) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<wgt_indx_bits, 1LL<<wgt_indx_bits, -1, 
                           wgt_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("wgt backend = %d, norm strat = %d, wgt indx bits = %d:\n", 
            wgt_backend, norm_strat, wgt_indx_bits);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_wgt_storage_types; backend++) {
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

int main()
{
    return runtests();
}
