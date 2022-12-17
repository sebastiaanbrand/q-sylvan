#include <inttypes.h>
#include <stdio.h>
#include <time.h>

#include "test_assert.h"
#include "qsylvan.h"

#include "grover.h"
#include "grover_cnf.h"
#include "shor.h"

bool VERBOSE = true;
double TOLERANCE = 1e-14;


int test_grover()
{
    BDDVAR nbits;
    AMP a;
    QMDD grov;
    double prob;

    // 2 qubit test (+1 ancilla)
    nbits = 2;
    bool flag2[] = {0,1}; // "flagged" entry
    grov = qmdd_grover(nbits, flag2);
    test_assert(qmdd_is_unitvector(grov, nbits+1));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag2[0] && x[1] == flag2[1])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(a == AADD_ZERO);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(fabs(prob - 1) < TOLERANCE);

    if(VERBOSE) printf("qmdd %2d-qubit Grover:        ok (Pr(flag) = %lf)\n", nbits+1, prob);


    // 3 qubit test (+1 ancilla)
    nbits = 3;
    bool flag3[] = {1,1,0}; // "flagged" entry
    grov = qmdd_grover(nbits, flag3);
    test_assert(qmdd_is_unitvector(grov, nbits+1));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag3[0] && x[1] == flag3[1] && x[2] == flag3[2])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(qmdd_amp_to_prob(a) < 0.004);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(prob > 0.94);

    if(VERBOSE) printf("qmdd %2d-qubit Grover:        ok (Pr(flag) = %lf)\n", nbits+1, prob);

    // 10 qubit test (+1 ancilla)
    nbits = 10;
    bool flag10[nbits];
    srand(time(NULL));
    for (BDDVAR i = 0; i < nbits; i++) flag10[i] = (bool)(rand() % 2);
    grov = qmdd_grover(nbits, flag10);
    test_assert(qmdd_is_close_to_unitvector(grov, nbits+1, 1e-12));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag10[0] && x[1] == flag10[1] && 
            x[2] == flag10[2] && x[3] == flag10[3] &&
            x[4] == flag10[4] && x[5] == flag10[5] &&
            x[6] == flag10[6] && x[7] == flag10[7] &&
            x[8] == flag10[8] && x[9] == flag10[9])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(qmdd_amp_to_prob(a) < 0.001);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(prob > 0.99);

    if(VERBOSE) printf("qmdd %2d-qubit Grover:        ok (Pr(flag) = %lf)\n", nbits+1, prob);
    return 0;
}

int test_grover_matrix()
{
    aadd_set_auto_gc_wgt_table(false); // no auto gc of ctable yet for mult operations

    BDDVAR nbits;
    AMP a;
    QMDD grov;
    double prob;

    // 2 qubit test (+1 ancilla)
    nbits = 2;
    bool flag2[] = {0,1}; // "flagged" entry
    grov = qmdd_grover_matrix(nbits, flag2);
    test_assert(qmdd_is_unitvector(grov, nbits+1));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag2[0] && x[1] == flag2[1])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(a == AADD_ZERO);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(fabs(prob - 1) < TOLERANCE);

    if(VERBOSE) printf("matrix qmdd %2d-qubit Grover: ok (Pr(flag) = %lf)\n", nbits+1, prob);


    // 3 qubit test (+1 ancilla)
    nbits = 3;
    bool flag3[] = {1,1,0}; // "flagged" entry
    grov = qmdd_grover_matrix(nbits, flag3);
    test_assert(qmdd_is_unitvector(grov, nbits+1));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag3[0] && x[1] == flag3[1] && x[2] == flag3[2])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(qmdd_amp_to_prob(a) < 0.004);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(prob > 0.94);

    if(VERBOSE) printf("matrix qmdd %2d-qubit Grover: ok (Pr(flag) = %lf)\n", nbits+1, prob);


    // 10 qubit test (+1 ancilla)
    nbits = 10;
    bool flag10[nbits];
    srand(time(NULL));
    for (BDDVAR i = 0; i < nbits; i++) flag10[i] = (bool)(rand() % 2);
    grov = qmdd_grover_matrix(nbits, flag10);
    test_assert(qmdd_is_close_to_unitvector(grov, nbits+1, TOLERANCE*100));

    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nbits+1)); k++) {
        bool *x = int_to_bitarray(k, nbits+1, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == flag10[0] && x[1] == flag10[1] && 
            x[2] == flag10[2] && x[3] == flag10[3] &&
            x[4] == flag10[4] && x[5] == flag10[5] &&
            x[6] == flag10[6] && x[7] == flag10[7] &&
            x[8] == flag10[8] && x[9] == flag10[9])
            prob += qmdd_amp_to_prob(a);
        else
            test_assert(qmdd_amp_to_prob(a) < 0.001);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(prob > 0.99);

    if(VERBOSE) printf("matrix qmdd %2d-qubit Grover: ok (Pr(flag) = %lf)\n", nbits+1, prob);
    return 0;

    aadd_set_auto_gc_wgt_table(true);
    return 0;
}

int test_grover_cnf()
{
    BDDVAR nqubits;
    BDDVAR k;
    BDDVAR clauses;
    BDDVAR answers;
    AMP a;
    QMDD grov;
    double prob;


    // 3-SAT, 3 qubit test
    nqubits = 3;
    // "flagged" cnf (~x1 V ~x2 V ~x3)^(~x1 V ~x2 V x3)^(~x1 V x2 V ~x3)^(x1 V ~x2 V ~x3)^(x1 V ~x2 V x3)^(x1 V x2 V ~x3)^(x1 V x2 V x3)
    int cnf3[] = {-1,-2,-3, -1,-2,3, -1,2,-3, 1,-2,-3, -1,2,3, 1,-2,3, 1,2,-3}; // In the order c1(x3,x2,x1),c2(x3,x2,x1), etc.
    bool ans3[] = {1,1,1};
    k = 3;
    clauses = 7;
    answers = 1;
    grov = qmdd_grover_cnf(nqubits, cnf3, k, clauses, answers);
    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nqubits+1+clauses)); k++) {
        bool *x = int_to_bitarray(k, nqubits+1+clauses, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == ans3[0] && x[1] == ans3[1])
            prob += qmdd_amp_to_prob(a);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(qmdd_is_unitvector(grov, nqubits+1+clauses));

    if(VERBOSE) printf("qmdd %2d-qubit 3-SAT Grover:  ok (Pr(flag) = %lf)\n", nqubits, prob);


    // 3-SAT, 4 qubit test
    nqubits = 4;
    // "flagged" cnf (~x1 V ~x2 V ~x3)^(~x1 V ~x2 V x3)^(~x1 V x2 V ~x3)^(x1 V ~x2 V ~x3)^(x1 V ~x2 V x3)^(x1 V x2 V ~x3)^(x1 V x2 V x3)^(~X1 V ~X2 V X4)
    int cnf4[] = {-1,-2,-3, -1,-2,3, -1,2,-3, 1,-2,-3, -1,2,3, 1,-2,3, 1,2,-3, 1,2,-4}; // In the order c1(x3,x2,x1),c2(x3,x2,x1), etc.
    bool ans4[] = {1,1,1,1};
    k = 3;
    clauses = 8;
    answers = 1;
    grov = qmdd_grover_cnf(nqubits, cnf4, k, clauses, answers);
    // test probabilities
    prob = 0;
    for (int k = 0; k < (1<<(nqubits+1+clauses)); k++) {
        bool *x = int_to_bitarray(k, nqubits+1+clauses, true);
        a = aadd_getvalue(grov, x);
        if (x[0] == ans4[0] && x[1] == ans4[1] && x[2] == ans4[2] && x[3] == ans4[3])
            prob += qmdd_amp_to_prob(a);
        free(x); // int_to_bitarray mallocs
    }
    test_assert(qmdd_is_unitvector(grov, nqubits+1+clauses));

    if(VERBOSE) printf("qmdd %2d-qubit 3-SAT Grover:  ok (Pr(flag) = %lf)\n", nqubits, prob);

    return 0;
}

int test_shor()
{   
    qmdd_shor_set_testing_mode(true); // internal sanity checks in shor implementation
    QMDD q, qref;
    bool x3[] = {0,0,0};
    bool as[] = {1,0,1};
    AMP a;

    // <Test qmdd_phi_add>
    // Test inversion
    // (no controls)
    BDDVAR nc = AADD_INVALID_VAR;
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_gate(q, GATEID_H, 0);
    q = qmdd_gate(q, GATEID_H, 1);
    q = qmdd_gate(q, GATEID_H, 2);
    qref = q;
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);
    q = qmdd_phi_add_inv(q, 0, 2, nc, nc, as);
    test_assert(aadd_equivalent(q, qref, 3, false, false));
    test_assert(aadd_equivalent(q, qref, 3, true, false));
    test_assert(q == qref);

    // Test addition in Fourier space (be mindful about endianness here!)
    // 2 + 1 (no carry)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 0; as[2] = 0;          // a = 100 = 1 (LSB first)
    q = qmdd_create_basis_state(3, x3);        // create state |x>
    q = qmdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);     // addition in Fourier space gives |phi(x+a)>
    q = qmdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |3> = |011> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[0]=1; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    // 2 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qmdd_create_basis_state(3, x3);        // create state |x>
    q = qmdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);     // addition in Fourier space gives |phi(x+a)>
    q = qmdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |4> = |100> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[0]=1; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    // 2 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qmdd_create_basis_state(3, x3);        // create state |x>
    q = qmdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);     // addition in Fourier space gives |phi(x+a)>
    q = qmdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    // 3 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qmdd_create_basis_state(3, x3);        // create state |x>
    q = qmdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);     // addition in Fourier space gives |phi(x+a)>
    q = qmdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    // 3 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qmdd_create_basis_state(3, x3);        // create state |x>
    q = qmdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qmdd_phi_add(q, 0, 2, nc, nc, as);     // addition in Fourier space gives |phi(x+a)>
    q = qmdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |6> = |110> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[0]=1; x3[1]=1; x3[2]=1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    // </Test qmdd_phi_add>

    if(VERBOSE) printf("qmdd fourier addition:       ok\n");


    // test Shor
    srand(time(NULL));
    uint64_t N, factor, counter, nqubits;

    // 15 = 3 x 5 (11 qubits)
    N = 15;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = shor_run(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qmdd %" PRIu64 "-qubit Shor:          ok (found factor %" PRIu64 " of %" PRIu64 " with %" PRIu64 " tries)\n", nqubits, factor, N, counter);

    // 35 = 5 x 7 (15 qubits)
    N = 35;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = shor_run(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qmdd %" PRIu64 "-qubit Shor:          ok (found factor %" PRIu64 " of %" PRIu64 " with %" PRIu64 " tries)\n", nqubits, factor, N, counter);


    return 0;
}

int runtests()
{
    if (test_grover()) return 1;
    if (test_grover_matrix()) return 1;
    if (test_grover_cnf()) return 1;
    if (test_shor()) return 1;

    return 0;
}

int main()
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s)\n", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<16, TOLERANCE, COMP_HASHMAP, NORM_LARGEST);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    int res = runtests();

    sylvan_quit();

    return res;
}
