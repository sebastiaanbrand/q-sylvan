#include "supremacy.h"

static uint32_t
random_sqrtXY()
{
    return ( rand() % 2 ) ? GATEID_sqrtX : GATEID_sqrtY;
}

QDD
supremacy_5_1_circuit(uint32_t depth)
{
    LACE_ME;

    // Start with |00000> state
    BDDVAR n_qubits = 5;
    QDD qdd = qdd_create_all_zero_state(n_qubits);

    // H on all qubits
    for (BDDVAR k=0; k < n_qubits; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // First 4 cycles the single qubit gates we apply are all T gates
    if (depth > 0) {
        qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);
        qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);
    }
    if (depth > 1) {
        qdd = qdd_cgate(qdd, GATEID_Z, 1, 2);
        qdd = qdd_gate(qdd, GATEID_T, 0);
        qdd = qdd_gate(qdd, GATEID_T, 3);
        qdd = qdd_gate(qdd, GATEID_T, 4);
    }

    if (depth > 2) {
        qdd = qdd_cgate(qdd, GATEID_Z, 2, 3);
        qdd = qdd_gate(qdd, GATEID_T, 1);
    }
    if (depth > 3) {
        qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);
        qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);
        qdd = qdd_gate(qdd, GATEID_T, 2);
    }

    // Following cycles the single qubit gates are random from {sqrt(X), sqrt(Y)}
    for (uint32_t d = 4; d < depth; d++) {
        switch (d % 3) {
        case 0:
            qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);    // CZ(0,1)
            qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);    // CZ(3,4)
            qdd = qdd_gate(qdd, random_sqrtXY(), 2); // random on q2
            break;
        case 1:
            qdd = qdd_cgate(qdd, GATEID_Z, 1, 2);    // CZ(1,2)
            qdd = qdd_gate(qdd, random_sqrtXY(), 0); // random on q0
            qdd = qdd_gate(qdd, random_sqrtXY(), 3); // random on q3
            qdd = qdd_gate(qdd, random_sqrtXY(), 4); // random on q4
            break;
        case 2:
            qdd = qdd_cgate(qdd, GATEID_Z, 2, 3);    // CZ(2,3)
            qdd = qdd_gate(qdd, random_sqrtXY(), 0); // random on q0
            qdd = qdd_gate(qdd, random_sqrtXY(), 1); // random on q1
            qdd = qdd_gate(qdd, random_sqrtXY(), 4); // random on q4
            break;
        default:
            break;
        }
    }
    return qdd;
}

QDD
supremacy_5_4_circuit(uint32_t depth)
{
    // 5x4 "grid" from [Characterizing Quantum Supremacy in Near-Term Devices]
    // (maybe a lot of hard-coded stuff but oh well...)
    LACE_ME;

    // Start with |00...0> state
    BDDVAR n_qubits = 20;
    QDD qdd = qdd_create_all_zero_state(n_qubits);

    // H on all qubits
    for (BDDVAR k=0; k < n_qubits; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // qubits which are not involved in a CZ at cycle (d mod 8) 
    // AND had a CZ applied to them in the previous cycle.
    BDDVAR qubits0[6]  = {1,8,10,14,17,19};
    BDDVAR qubits1[8]  = {2,3,5,6,12,13,15,16};
    BDDVAR qubits2[6]  = {0,1,7,10,17,18};
    BDDVAR qubits3[4]  = {6,8,11,13};
    BDDVAR qubits4[4]  = {5,9,10,12};
    BDDVAR qubits5[8]  = {3,4,6,7,13,14,16,17};
    BDDVAR qubits6[4]  = {1,8,12,19};
    BDDVAR qubits7[10] = {0,2,4,5,7,9,11,13,16,18};

    // First 7 cycles the single qubit gates we apply are all T gates
    // (not sure if we also should do random {sqrt(X), sqrt(Y)} on the
    // remaining qubits, the paper is not super clear)
    if (depth > 0) {
        qdd = qdd_cgate(qdd, GATEID_Z,  2,  3);    // CZ(2,3)
        qdd = qdd_cgate(qdd, GATEID_Z,  5,  6);    // CZ(5,6)
        qdd = qdd_cgate(qdd, GATEID_Z, 12, 13);    // CZ(12,13)
        qdd = qdd_cgate(qdd, GATEID_Z, 15, 16);    // CZ(15,16)
    }
    if (depth > 1) {
        qdd = qdd_cgate(qdd, GATEID_Z,  0,  1);    // CZ(0,1)
        qdd = qdd_cgate(qdd, GATEID_Z,  7,  8);    // CZ(7,8)
        qdd = qdd_cgate(qdd, GATEID_Z, 10, 11);    // CZ(10,11)
        qdd = qdd_cgate(qdd, GATEID_Z, 17, 18);    // CZ(17,18)
        BDDVAR qubits[8] = {2,3,5,6,12,13,15,16};
        for (int i = 0; i < 8; i++)
            qdd = qdd_gate(qdd, GATEID_T, qubits[i]);
    }
    if (depth > 2) {
        qdd = qdd_cgate(qdd, GATEID_Z,  6, 11);    // CZ(6,11)
        qdd = qdd_cgate(qdd, GATEID_Z,  8, 13);    // CZ(8,13)
        BDDVAR qubits[8] = {0,1,7,8,10,11,17,18};
        for (int i = 0; i < 8; i++)
            qdd = qdd_gate(qdd, GATEID_T, qubits[i]);
    }
    if (depth > 3) {
        qdd = qdd_cgate(qdd, GATEID_Z,  5, 10);    // CZ(5,10)
        qdd = qdd_cgate(qdd, GATEID_Z,  7, 12);    // CZ(7,12)
        qdd = qdd_cgate(qdd, GATEID_Z,  9, 14);    // CZ(9,14)
    }
    if (depth > 4) {
        qdd = qdd_cgate(qdd, GATEID_Z,  3,  4);    // CZ(3,4)
        qdd = qdd_cgate(qdd, GATEID_Z,  6,  7);    // CZ(6,7)
        qdd = qdd_cgate(qdd, GATEID_Z, 13, 14);    // CZ(13,14)
        qdd = qdd_cgate(qdd, GATEID_Z, 16, 17);    // CZ(16,17)
        qdd = qdd_gate(qdd, GATEID_T,  9);
        qdd = qdd_gate(qdd, GATEID_T, 14);
    }
    if (depth > 5) {
        qdd = qdd_cgate(qdd, GATEID_Z,  1,  2);    // CZ(1,2)
        qdd = qdd_cgate(qdd, GATEID_Z,  8,  9);    // CZ(8,9)
        qdd = qdd_cgate(qdd, GATEID_Z, 11, 12);    // CZ(11,12)
        qdd = qdd_cgate(qdd, GATEID_Z, 18, 19);    // CZ(18,19)
        qdd = qdd_gate(qdd, GATEID_T, 4);
    }
    if (depth > 6) {
        qdd = qdd_cgate(qdd, GATEID_Z,  0,  5);    // CZ(0,5)
        qdd = qdd_cgate(qdd, GATEID_Z,  2,  7);    // CZ(2,7)
        qdd = qdd_cgate(qdd, GATEID_Z,  4,  9);    // CZ(4,9)
        qdd = qdd_cgate(qdd, GATEID_Z, 11, 16);    // CZ(11,16)
        qdd = qdd_cgate(qdd, GATEID_Z, 13, 18);    // CZ(13,18)
        qdd = qdd_gate(qdd, GATEID_T, 19);
    }
    // Following cycles the single qubit gates are random from {sqrt(X), sqrt(Y)}
    for (uint32_t d = 7; d <= depth; d++) {
        switch (d % 8) {
        case 0:
            qdd = qdd_cgate(qdd, GATEID_Z,  2,  3);    // CZ(2,3)
            qdd = qdd_cgate(qdd, GATEID_Z,  5,  6);    // CZ(5,6)
            qdd = qdd_cgate(qdd, GATEID_Z, 12, 13);    // CZ(12,13)
            qdd = qdd_cgate(qdd, GATEID_Z, 15, 16);    // CZ(15,16)
            for (int i = 0; i < 6; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits0[i]);
            break;
        case 1:
            qdd = qdd_cgate(qdd, GATEID_Z,  0,  1);    // CZ(0,1)
            qdd = qdd_cgate(qdd, GATEID_Z,  7,  8);    // CZ(7,8)
            qdd = qdd_cgate(qdd, GATEID_Z, 10, 11);    // CZ(10,11)
            qdd = qdd_cgate(qdd, GATEID_Z, 17, 18);    // CZ(17,18)
            for (int i = 0; i < 8; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits1[i]);
            break;
        case 2:
            qdd = qdd_cgate(qdd, GATEID_Z,  6, 11);    // CZ(6,11)
            qdd = qdd_cgate(qdd, GATEID_Z,  8, 13);    // CZ(8,13)
            for (int i = 0; i < 6; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits2[i]);
            break;
        case 3:
            qdd = qdd_cgate(qdd, GATEID_Z,  5, 10);    // CZ(5,10)
            qdd = qdd_cgate(qdd, GATEID_Z,  7, 12);    // CZ(7,12)
            qdd = qdd_cgate(qdd, GATEID_Z,  9, 14);    // CZ(9,14)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits3[i]);
            break;
        case 4:
            qdd = qdd_cgate(qdd, GATEID_Z,  3,  4);    // CZ(3,4)
            qdd = qdd_cgate(qdd, GATEID_Z,  6,  7);    // CZ(6,7)
            qdd = qdd_cgate(qdd, GATEID_Z, 13, 14);    // CZ(13,14)
            qdd = qdd_cgate(qdd, GATEID_Z, 16, 17);    // CZ(16,17)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits4[i]);
            break;
        case 5:
            qdd = qdd_cgate(qdd, GATEID_Z,  1,  2);    // CZ(1,2)
            qdd = qdd_cgate(qdd, GATEID_Z,  8,  9);    // CZ(8,9)
            qdd = qdd_cgate(qdd, GATEID_Z, 11, 12);    // CZ(11,12)
            qdd = qdd_cgate(qdd, GATEID_Z, 18, 19);    // CZ(18,19)
            for (int i = 0; i < 8; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits5[i]);
            break;
        case 6:
            qdd = qdd_cgate(qdd, GATEID_Z,  0,  5);    // CZ(0,5)
            qdd = qdd_cgate(qdd, GATEID_Z,  2,  7);    // CZ(2,7)
            qdd = qdd_cgate(qdd, GATEID_Z,  4,  9);    // CZ(4,9)
            qdd = qdd_cgate(qdd, GATEID_Z, 11, 16);    // CZ(11,16)
            qdd = qdd_cgate(qdd, GATEID_Z, 13, 18);    // CZ(13,18)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits6[i]);
            break;
        case 7:
            qdd = qdd_cgate(qdd, GATEID_Z,  1,  6);    // CZ(1,6)
            qdd = qdd_cgate(qdd, GATEID_Z,  3,  8);    // CZ(3,8)
            qdd = qdd_cgate(qdd, GATEID_Z, 10, 15);    // CZ(10,15)
            qdd = qdd_cgate(qdd, GATEID_Z, 12, 17);    // CZ(12,17)
            qdd = qdd_cgate(qdd, GATEID_Z, 14, 19);    // CZ(14,19)
            for (int i = 0; i < 10; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits7[i]);
            break;
        default:
            break;
        }
    }
    return qdd;
}
