#include <stdio.h>

#include "qsylvan.h"
#include "test_assert.h"
#include "../examples/grover.h"

int test_grover_gc()
{
    int qubits = 8;
    bool *flag = qmdd_grover_ones_flag(qubits+1); // + 1 ancilla

    // run Grover
    QMDD qmdd = qmdd_grover(qubits, flag);

    // Sanity checks on final state
    flag[qubits] = 0; AADD_WGT amp0 = aadd_getvalue(qmdd, flag);
    flag[qubits] = 1; AADD_WGT amp1 = aadd_getvalue(qmdd, flag);
    double flag_prob = qmdd_amp_to_prob(amp0) + qmdd_amp_to_prob(amp1);
    free(flag);

    test_assert(flag_prob > 0.9 && flag_prob < 1.0+1e-6);

    return 0;
}


int run_qmdd_tests()
{
    // Test gc by running some circuits for which gc triggers
    if (test_grover_gc()) return 1;

    return 0;
}


int test_with(int amps_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Initialize Q-Sylvan with tolerance 0 (this creates larger QMDDs such that
    // garbage collection of the edge weight table is triggered earlier)
    double tol = 0;
    uint64_t wgt_tab_size = 1LL<<15;
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, tol, amps_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("amps backend = %d, norm strategy = %d:\n", amps_backend, norm_strat);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_exit();

    return res;
}

int runtests()
{
    int backend = COMP_HASHMAP;
    for (int norm_strat = 0; norm_strat < 2; norm_strat++) {
        if (test_with(backend, norm_strat)) return 1;
    }
    return 0;
}

int main()
{
    return runtests();
}
