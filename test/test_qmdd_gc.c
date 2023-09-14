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


int test_with(int wgt_backend, int norm_strat, int wgt_indx_bits) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Initialize Q-Sylvan with tolerance 0 (this creates larger QMDDs such that
    // garbage collection of the edge weight table is triggered earlier)
    double tol = 0;
    uint64_t wgt_tab_size = 1LL<<wgt_indx_bits;
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, wgt_tab_size, tol, wgt_backend, norm_strat);
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
    int backend = COMP_HASHMAP;
    for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
        if (test_with(backend, norm_strat, 15)) return 1;
        if (backend == COMP_HASHMAP) {
            // test with edge wgt index > 23 bits
            if (test_with(backend, norm_strat, 24)) return 1;
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
