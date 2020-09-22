#include <stdio.h>
#include <time.h>

#include "test_assert.h"
#include "sylvan.h"
#include "sylvan_qdd_complex.h"

#include "grover.h"

bool VERBOSE = true;

int test_grover()
{
    BDDVAR nqubits;
    AMP a;
    QDD grov;
    double prob;

    // 2 qubit test
    nqubits = 2;
    bool x2[] = {0,1}; // "flagged" entry
    grov = qdd_grover(nqubits, x2);
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(grov, x2); test_assert(a == C_ZERO);
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(grov, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(grov, x2); test_assert(a == C_ONE);  prob = comp_to_prob(comp_value(a));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(grov, x2); test_assert(a == C_ZERO);
    test_assert(qdd_is_unitvector(grov, nqubits));

    if(VERBOSE) printf("qdd %2d-qubit Grover:      ok (Pr(flag) = %lf)\n", nqubits, prob);


    // 3 qubit test
    nqubits = 3;
    bool x3[] = {1,1,0}; // "flagged" entry
    grov = qdd_grover(3, x3);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) > 0.94);  prob = comp_to_prob(comp_value(a));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(grov, x3); test_assert(comp_to_prob(comp_value(a)) < 0.008);
    test_assert(qdd_is_unitvector(grov, 3));

    if(VERBOSE) printf("qdd %2d-qubit Grover:      ok (Pr(flag) = %lf)\n", nqubits, prob);


    // 10 qubit test (random flag)
    nqubits = 10;
    bool x10[nqubits];
    srand(time(NULL));
    for (BDDVAR i = 0; i < nqubits; i++) x10[i] = (bool)(rand() % 2);
    grov = qdd_grover(nqubits, x10);
    test_assert(qdd_is_close_to_unitvector(grov, nqubits, TOLERANCE*100));
    prob = comp_to_prob(comp_value(qdd_get_amplitude(grov, x10)));

    if(VERBOSE) printf("qdd %2d-qubit Grover:      ok (Pr(flag) = %lf)\n", nqubits, prob);
    return 0;
}

int test_shor()
{   
    QDD q, qref;
    bool x3[] = {0,0,0};
    bool as[] = {1,0,1};
    AMP a;

    LACE_ME;

    // <Test qdd_phi_add>
    // Test inversion
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0);
    q = qdd_gate(q, GATEID_H, 1);
    q = qdd_gate(q, GATEID_H, 2);
    qref = q;
    q = qdd_phi_add(q, 0, 2, as);
    q = qdd_phi_add_inv(q, 0, 2, as);
    test_assert(qdd_equivalent(q, qref, 3, false, false));
    test_assert(qdd_equivalent(q, qref, 3, true, false));
    test_assert(q == qref);

    // Test addition in Fourier space (be mindful about endianness here!)
    // 2 + 1 (no carry)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 0; as[2] = 0;          // a = 100 = 1 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |3> = |011> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 2 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |4> = |100> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 2 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 3 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 3 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |6> = |110> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    // </Test qdd_phi_add>

    if(VERBOSE) printf("qdd fourier addition:     ok\n");


    // test Shor
    srand(time(NULL));
    uint64_t N, factor, counter, nqubits;

    // 15 = 3 x 5 (11 qubits)
    N = 15;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = run_shor(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qdd %ld-qubit Shor:        ok (found factor %ld of %ld with %ld tries)\n", nqubits, factor, N, counter);

    // 35 = 5 x 7 (15 qubits)
    N = 35;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = run_shor(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qdd %ld-qubit Shor:        ok (found factor %ld of %ld with %ld tries)\n", nqubits, factor, N, counter);


    return 0;
}

int runtests()
{
    if (test_grover()) return 1;
    if (test_shor()) return 1;

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

    sylvan_quit();
    lace_exit();

    return res;
}