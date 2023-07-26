#include <qsylvan.h>
#include <stdlib.h> 

void ry_cz_ansatz(int nqubits, int depth)
{
    QMDD state = qmdd_create_all_zero_state(nqubits);
    aadd_protect(&state);

    for (int d = 0; d < depth; d++) {
        // Ry rotations
        for (int n = 0; n < nqubits; n++) {
            fl_t theta = (fl_t)rand() / (fl_t)RAND_MAX;
            state = qmdd_gate(state, GATEID_Ry(theta), n);
        }

        // entangling gates
        for (int n = 0; n < nqubits - 1; n++) {
            state = qmdd_cgate(state, GATEID_Z, n, n+1);
        }
    }

    // measure
    bool *outcome = malloc(sizeof(bool) * nqubits);
    double prob;
    qmdd_measure_all(state, nqubits, outcome, &prob);
    printf("measured state: |");
    for (int n = 0; n < nqubits; n++){
        printf("%d", outcome[n]);
    }
    printf("> (with prob %lf)\n", prob);

    free(outcome);
    aadd_unprotect(&state);
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

    srand(time(NULL));
    ry_cz_ansatz(10, 5);

    sylvan_quit();
    return 0;
}
