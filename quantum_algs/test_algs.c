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

int runtests()
{
    if (test_grover()) return 1;

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