#include <stdio.h>
#include <sys/time.h>

#include "sylvan.h"
#include "grover.h"
#include "supremacy.h"

/**
 * Obtain current wallclock time
 */
static double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

int bench_supremacy_5_1(uint32_t depth, uint32_t workers)
{
    // 5x1 "grid" from [Characterizing Quantum Supremacy in Near-Term Devices]
    printf("bench sup5_1, depth %d, %2d worker(s)\n", depth, workers);

    uint64_t node_count;
    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);
    LACE_ME;

    // Init Sylvan
    sylvan_set_limits(4LL<<30, 1, 6);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<23);

    srand(42);
    QDD res = supremacy_5_1_circuit(depth);
    node_count = qdd_countnodes(res);

    t_end = wctime();
    runtime = (t_end - t_start);

    printf("%lf sec", runtime);
    printf(", %ld nodes", node_count);
    printf("\n");

    // Cleanup
    sylvan_quit();
    lace_exit();

    return 0;
}

int
bench_supremacy_5_4(uint32_t depth, uint32_t workers)
{
    printf("bench sup5_4, depth %d, %2d worker(s)\n", depth, workers);

    uint64_t node_count;
    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_limits(4LL<<30, 1, 6);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<23);

    srand(66);

    QDD res = supremacy_5_4_circuit(depth);
    node_count = qdd_countnodes(res); 

    t_end = wctime();
    runtime = (t_end - t_start);

    printf("%lf sec", runtime);
    printf(", %ld nodes", node_count);
    printf("\n");

    // Cleanup
    sylvan_quit();
    lace_exit();

    return 0;
}


int bench_grover(int num_qubits, bool flag[], int workers, FILE *logfile)
{
    printf("bench grover, %d qubits, %2d worker(s), ", num_qubits, workers); 
    printf("flag = [");
    for (int i = 0; i < num_qubits; i++)
        printf("%d",flag[i]);
    printf("], ");
    fflush(stdout);

    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_limits(4LL<<30, 1, 6);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<18);

    QDD grov;
    uint64_t node_count;
    
    if (logfile != NULL) qdd_stats_start(logfile);

    grov = qdd_grover(num_qubits, flag);

    if (logfile != NULL) qdd_stats_finish();

    t_end = wctime();
    runtime = (t_end - t_start);

    node_count = qdd_countnodes(grov); // TODO: also get Pr(flag)
    printf("%ld nodes, %lf sec \n", node_count, runtime);

    // Cleanup
    sylvan_quit();
    lace_exit();
    return 0;
}

int main()
{
    int n = 21;
    bool flag[] = {1,1,1,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,0,1};
    bench_grover(n, flag, 1, NULL);


    bench_supremacy_5_1(100, 1);
    bench_supremacy_5_4(5, 1);

    return 0;
}