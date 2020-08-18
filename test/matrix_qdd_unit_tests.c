#include <stdio.h>
#include <time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_int.h"


int test_matrix_creation()
{
    QDD matrix = qdd_create_single_qubit_gate(3, 1, GATEID_Z);
    FILE *fp;
    fp = fopen("matrix_test.dot", "w");
    qdd_fprintdot(fp, matrix, false);

    return 0;
}

int runtests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    if (test_matrix_creation()) return 1;

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

    free_amplitude_table();
    sylvan_quit();
    lace_exit();

    return res;
}
