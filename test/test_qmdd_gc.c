#include <stdio.h>

#include "qsylvan.h"
#include "test_assert.h"
#include "../examples/grover.h"

uint64_t min_wgt_tablesize = 1LL<<14;
uint64_t max_wgt_tablesize = 1LL<<17;


int test_grover_gc()
{
    int qubits = 8;
    bool *flag = qmdd_grover_ones_flag(qubits+1); // + 1 ancilla

    // run Grover
    QMDD qmdd = qmdd_grover(qubits, flag);

    // Sanity checks on final state
    flag[qubits] = 0; EVBDD_WGT amp0 = evbdd_getvalue(qmdd, flag);
    flag[qubits] = 1; EVBDD_WGT amp1 = evbdd_getvalue(qmdd, flag);
    double flag_prob = qmdd_amp_to_prob(amp0) + qmdd_amp_to_prob(amp1);
    free(flag);

    printf("flag prob = %lf\n", flag_prob);
    test_assert(flag_prob > 0.9 && flag_prob < 1.0+1e-6);

    return 0;
}


int run_qmdd_tests()
{
    // Test gc by running some circuits for which gc triggers
    if (test_grover_gc()) return 1;

    return 0;
}


int test_table_size_increase() 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);

    // Initialize Q-Sylvan with tolerance 0 (this creates larger QMDDs such that
    // garbage collection of the edge weight table is triggered earlier)
    double tol = 0;
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(min_wgt_tablesize, max_wgt_tablesize, -1, COMP_HASHMAP, NORM_MAX);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    // check that wgt_tablesize doubles after gc, but not beyond max_wgt_tablesize
    uint64_t wgt_tablesize = min_wgt_tablesize;
    for (int i = 0; i < 10; i++) {
        test_assert(sylvan_get_edge_weight_table_size() == wgt_tablesize);
        evbdd_gc_wgt_table();
        wgt_tablesize = 2*wgt_tablesize;
        if (wgt_tablesize > max_wgt_tablesize) {
            wgt_tablesize = max_wgt_tablesize;
        }
    }

    sylvan_quit();
    lace_stop();
    return 0;
}


int test_custom_gate_gc_protection()
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);

    // Initialize Q-Sylvan with tolerance 0 (this creates larger QMDDs such that
    // garbage collection of the edge weight table is triggered earlier)
    double tol = 1e-14;
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(min_wgt_tablesize, max_wgt_tablesize, -1, COMP_HASHMAP, NORM_MAX);
    qmdd_set_testing_mode(true); // turn on internal sanity tests


    QMDD qRef, qTest, qInit;
    BDDVAR nqubits, t;
    int tmp_gateid;
    double pi = 2.0 * flt_acos(0.0);

    nqubits = 5, t = 3;
    qInit = qmdd_create_all_zero_state(nqubits);

    // apply the same custom gate twice (w/o gc)
    qRef  = qmdd_gate(qInit, GATEID_U(pi/2.0, -pi/2.0, pi/4.0), t);
    qTest = qmdd_gate(qInit, GATEID_U(pi/2.0, -pi/2.0, pi/4.0), t);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, true));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // apply the same custom gate twice, but with gc of node table between them
    evbdd_protect(&qRef);
    evbdd_protect(&qTest);
    qRef  = qmdd_gate(qInit, GATEID_U(pi/2.0, -pi/2.0, pi/4.0), t);
    sylvan_gc();
    qTest = qmdd_gate(qInit, GATEID_U(pi/2.0, -pi/2.0, pi/4.0), t);
    evbdd_unprotect(&qRef);
    evbdd_unprotect(&qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, true));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // trigger gc of node table after defining temp gate, but before applying it
    evbdd_protect(&qRef);
    evbdd_protect(&qTest);
    tmp_gateid = GATEID_U(pi/2.0, -pi/2.0, pi/4.0);
    qRef  = qmdd_gate(qInit, tmp_gateid, t);
    sylvan_gc();
    qTest = qmdd_gate(qInit, tmp_gateid, t);
    evbdd_unprotect(&qRef);
    evbdd_unprotect(&qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, true));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);

    // trigger gc of edge wgt table after defining temp gate, but before applying it
    evbdd_protect(&qRef);
    evbdd_protect(&qTest);
    tmp_gateid = GATEID_U(pi/2.0, -pi/2.0, pi/4.0);
    qRef  = qmdd_gate(qInit, tmp_gateid, t);
    evbdd_gc_wgt_table();
    qTest = qmdd_gate(qInit, tmp_gateid, t);
    evbdd_unprotect(&qRef);
    evbdd_unprotect(&qTest);
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, false, true));
    test_assert(evbdd_equivalent(qRef, qTest, nqubits, true, false));
    test_assert(qTest == qRef);


    sylvan_quit();
    lace_stop();
    return 0;
}


int test_with(int wgt_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Initialize Q-Sylvan with tolerance 0 (this creates larger QMDDs such that
    // garbage collection of the edge weight table is triggered earlier)
    double tol = 0;
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(min_wgt_tablesize, max_wgt_tablesize, tol, wgt_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("wgt backend = %d, norm strat = %d:\n", wgt_backend, norm_strat);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();
    return res;
}

int runtests()
{
    int backend = COMP_HASHMAP;
    for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
        if (test_with(backend, norm_strat)) return 1;
    }
    if (test_table_size_increase()) return 1;
    if (test_custom_gate_gc_protection()) return 1;
    return 0;
}

int main()
{
    return runtests();
}
