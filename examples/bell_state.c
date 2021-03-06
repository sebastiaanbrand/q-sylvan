#include <qsylvan.h>

void sample_bell_state()
{
    LACE_ME; // required to run Lace functions
    srand(time(NULL));

    // Create |Phi^+>
    int nqubits = 2;
    QDD state = qdd_create_all_zero_state(nqubits);
    state = qdd_gate(state, GATEID_H, 0);
    state = qdd_cgate(state, GATEID_X, 0, 1);
    state = qdd_gate(state, GATEID_X, 0);

    // Measure state
    bool outcome[] = {0, 0};
    double prob;
    qdd_measure_all(state, nqubits, outcome, &prob);

    printf("measured state: |%d%d> (with prob %.3lf)\n", outcome[0], outcome[1], prob);
}

int main() 
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd_defaults(1LL<<20);

    sample_bell_state();

    sylvan_quit();
    lace_exit();
    return 0;
}


