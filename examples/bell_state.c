#include <qsylvan.h>

void sample_bell_state()
{
    srand(time(NULL));

    // Create |Phi^+>
    int nqubits = 2;
    QMDD state = qmdd_create_all_zero_state(nqubits);
    state = qmdd_gate(state, GATEID_H, 0);
    state = qmdd_cgate(state, GATEID_X, 0, 1);
    state = qmdd_gate(state, GATEID_X, 0);

    // Measure state
    bool outcome[] = {0, 0};
    double prob;
    qmdd_measure_all(state, nqubits, outcome, &prob);

    printf("measured state: |%d%d> (with prob %.3lf)\n", outcome[0], outcome[1], prob);
}

int main() 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_defaults(1LL<<20);

    sample_bell_state();

    sylvan_quit();
    lace_stop();
    return 0;
}


