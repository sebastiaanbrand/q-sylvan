#include "random_circuit.h"

// don't use for heap / malloced arrays
#define len(x) (sizeof(x) / sizeof(x[0]))

void
random_qubit(BDDVAR nqubits, BDDVAR *t)
{
    *t = rand() % nqubits;
}

void
random_control_target(BDDVAR nqubits, BDDVAR *c, BDDVAR *t)
{
    *c = rand() % nqubits;
    do {
        *t = rand() % nqubits;
    }
    while (*t == *c);

    if (*c > *t) {
        BDDVAR temp;
        temp = *t;
        *t = *c;
        *c = temp;
    }
}

void
random_pauli(uint32_t *gateid)
{
    uint32_t pauli[] = {GATEID_X, GATEID_Y, GATEID_Z};
    uint32_t g = ( rand() % len(pauli) );
    *gateid = pauli[g];
}

void
random_cliff3(uint32_t *gateid)
{
    uint32_t cliff[] = {GATEID_X, GATEID_H, GATEID_S};
    uint32_t g = ( rand() % len(cliff) );
    *gateid = cliff[g];
}

void
random_cliff7(uint32_t *gateid)
{
    uint32_t cliff[] = {GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_sqrtX, GATEID_sqrtY, GATEID_S};
    uint32_t g = ( rand() % len(cliff) );
    *gateid = cliff[g];
}

void
random_univ(uint32_t *gateid)
{
    uint32_t univ[] = {GATEID_X, GATEID_Y, GATEID_Z, GATEID_H, GATEID_T};
    uint32_t g = ( rand() % len(univ) );
    *gateid = univ[g];
}


QMDD
qmdd_run_random_circuit(BDDVAR nqubits, uint64_t ngates, double cgate_ratio, uint64_t rseed)
{
    srand(rseed);

    LACE_ME;
    
    BDDVAR c, t;
    double r;
    uint32_t U;
    bool cgate;
    QMDD qmdd = qmdd_create_all_zero_state(nqubits);

    for (uint64_t g = 0; g < ngates; g++) {
        // choose 1 or 2 qubit gate
        r = ((double)rand() / (double)RAND_MAX);
        cgate = (r < cgate_ratio);
        if (cgate) {
            // choose target and control
            random_control_target(nqubits, &c, &t);
            qmdd = qmdd_cgate(qmdd, GATEID_Z, c, t);
        }
        else {
            // choose target qubit
            random_qubit(nqubits, &t);
            random_univ(&U); // TODO: parameterize gate-set
            qmdd = qmdd_gate(qmdd, U, t);
        }
    }

    return qmdd;
}

QMDD
qmdd_run_random_single_qubit_gates(BDDVAR nqubits, uint64_t ngates, uint64_t rseed)
{
    return qmdd_run_random_circuit(nqubits, ngates, 0.0, rseed);
}