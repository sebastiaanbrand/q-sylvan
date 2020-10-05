#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "sylvan.h"
#include "sylvan_qdd_complex.h"
#include "grover.h"
#include "supremacy.h"

#ifdef HAVE_PROFILER
#include <gperftools/profiler.h>
static char* profile_name = NULL; //"bench_qdd.prof";
#endif


static bool VERBOSE = false;

static size_t min_tablesize;
static size_t max_tablesize;
static size_t min_cachesize;
static size_t max_cachesize;
static size_t ctable_size;
static double ctable_tolerance;
static double ctable_gc_thres;

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

static void
write_parameters(FILE *file)
{
    fprintf(file, "{\n");
    fprintf(file, "  \"min_tablesize\": %ld,\n", min_tablesize);
    fprintf(file, "  \"max_tablesize\": %ld,\n", max_tablesize);
    fprintf(file, "  \"min_cachesize\": %ld,\n", min_cachesize);
    fprintf(file, "  \"max_cachesize\": %ld,\n", max_cachesize);
    fprintf(file, "  \"ctable_size\": %ld,\n", ctable_size);
    fprintf(file, "  \"ctable_tolerance\": %.5e,\n", ctable_tolerance);
    fprintf(file, "  \"ctable_gc_thres\": %lf,\n", ctable_gc_thres);
    fprintf(file, "  \"propagate_complex\": %d\n", propagate_complex);
    fprintf(file, "}\n");
    fclose(file);
}

int bench_25qubit_circuit(int workers)
{
    printf("qdd 25 qubit circuit: %2d worker(s), ", workers); 
    fflush(stdout);

    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);
    LACE_ME;

    // Init Sylvan
    sylvan_set_limits(4LL<<30, 1, 6);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<23, -1);

    QDD q;
    uint64_t node_count;
    bool x25[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


    // 25 qubit state
    q = qdd_create_basis_state(25, x25);

    // 500 gates
    q = qdd_cgate(q, GATEID_Z, 7, 21);	q = qdd_cgate(q, GATEID_Z, 4, 10);	q = qdd_cgate(q, GATEID_Z, 2, 12);	q = qdd_cgate(q, GATEID_Z, 4, 9);
    q = qdd_gate(q, GATEID_S, 4);    	q = qdd_cgate(q, GATEID_X, 10, 21);	q = qdd_gate(q, GATEID_Z, 16);    	q = qdd_gate(q, GATEID_H, 4);    
    q = qdd_gate(q, GATEID_S, 9);    	q = qdd_cgate(q, GATEID_Z, 4, 24);	q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_X, 2);    
    q = qdd_cgate(q, GATEID_X, 13, 23);	q = qdd_gate(q, GATEID_S, 6);    	q = qdd_cgate(q, GATEID_Z, 5, 22);	q = qdd_gate(q, GATEID_H, 24);    
    q = qdd_gate(q, GATEID_X, 17);    	q = qdd_cgate(q, GATEID_X, 14, 18);	q = qdd_gate(q, GATEID_H, 16);    	q = qdd_gate(q, GATEID_S, 23);    
    q = qdd_cgate(q, GATEID_Z, 13, 18);	q = qdd_cgate(q, GATEID_X, 0, 9);	q = qdd_gate(q, GATEID_S, 1);    	q = qdd_gate(q, GATEID_S, 14);    
    q = qdd_gate(q, GATEID_X, 21);    	q = qdd_cgate(q, GATEID_X, 8, 22);	q = qdd_cgate(q, GATEID_X, 10, 24);	q = qdd_cgate(q, GATEID_X, 15, 21);
    q = qdd_cgate(q, GATEID_X, 19, 23);	q = qdd_gate(q, GATEID_Z, 1);    	q = qdd_gate(q, GATEID_H, 3);    	q = qdd_cgate(q, GATEID_X, 1, 19);
    q = qdd_cgate(q, GATEID_X, 11, 17);	q = qdd_cgate(q, GATEID_Z, 14, 19);	q = qdd_cgate(q, GATEID_X, 6, 20);	q = qdd_gate(q, GATEID_X, 12);    
    q = qdd_cgate(q, GATEID_Z, 3, 22);	q = qdd_gate(q, GATEID_H, 7);    	q = qdd_cgate(q, GATEID_Z, 9, 14);	q = qdd_gate(q, GATEID_H, 13);    
    q = qdd_gate(q, GATEID_Z, 15);    	q = qdd_cgate(q, GATEID_X, 15, 22);	q = qdd_cgate(q, GATEID_X, 5, 20);	q = qdd_cgate(q, GATEID_X, 8, 14);
    q = qdd_cgate(q, GATEID_Z, 10, 17);	q = qdd_cgate(q, GATEID_X, 6, 18);	q = qdd_cgate(q, GATEID_X, 1, 22);	q = qdd_cgate(q, GATEID_Z, 21, 23);
    q = qdd_gate(q, GATEID_Z, 3);    	q = qdd_cgate(q, GATEID_X, 9, 18);	q = qdd_gate(q, GATEID_X, 18);    	q = qdd_cgate(q, GATEID_X, 2, 23);
    q = qdd_gate(q, GATEID_S, 17);    	q = qdd_cgate(q, GATEID_X, 12, 22);	q = qdd_gate(q, GATEID_H, 10);    	q = qdd_gate(q, GATEID_Z, 5);    
    q = qdd_gate(q, GATEID_S, 2);    	q = qdd_gate(q, GATEID_X, 19);    	q = qdd_cgate(q, GATEID_Z, 2, 22);	q = qdd_gate(q, GATEID_S, 15);    
    q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_X, 18);    	q = qdd_gate(q, GATEID_X, 16);    	q = qdd_cgate(q, GATEID_Z, 21, 23);
    q = qdd_cgate(q, GATEID_Z, 1, 13);	q = qdd_cgate(q, GATEID_X, 15, 18);	q = qdd_cgate(q, GATEID_X, 15, 22);	q = qdd_gate(q, GATEID_H, 14);    
    q = qdd_cgate(q, GATEID_X, 7, 16);	q = qdd_cgate(q, GATEID_Z, 8, 18);	q = qdd_cgate(q, GATEID_Z, 19, 20);	q = qdd_cgate(q, GATEID_Z, 14, 15);
    q = qdd_gate(q, GATEID_Z, 8);    	q = qdd_cgate(q, GATEID_X, 7, 14);	q = qdd_gate(q, GATEID_Z, 19);    	q = qdd_cgate(q, GATEID_Z, 17, 21);
    q = qdd_cgate(q, GATEID_X, 19, 20);	q = qdd_cgate(q, GATEID_X, 15, 20);	q = qdd_gate(q, GATEID_H, 2);    	q = qdd_cgate(q, GATEID_Z, 19, 23);
    q = qdd_cgate(q, GATEID_Z, 15, 20);	q = qdd_cgate(q, GATEID_Z, 2, 18);	q = qdd_gate(q, GATEID_X, 18);    	q = qdd_gate(q, GATEID_X, 6);    
    q = qdd_gate(q, GATEID_Z, 17);    	q = qdd_gate(q, GATEID_H, 24);    	q = qdd_cgate(q, GATEID_Z, 14, 16);	q = qdd_gate(q, GATEID_H, 9);    
    q = qdd_gate(q, GATEID_Z, 13);    	q = qdd_cgate(q, GATEID_Z, 7, 10);	q = qdd_gate(q, GATEID_Z, 2);    	q = qdd_cgate(q, GATEID_Z, 0, 15);
    q = qdd_cgate(q, GATEID_X, 7, 19);	q = qdd_cgate(q, GATEID_X, 21, 24);	q = qdd_gate(q, GATEID_Z, 7);    	q = qdd_gate(q, GATEID_Z, 1);    
    q = qdd_gate(q, GATEID_Z, 17);    	q = qdd_gate(q, GATEID_H, 3);    	q = qdd_cgate(q, GATEID_X, 20, 24);	q = qdd_cgate(q, GATEID_X, 4, 14);
    q = qdd_gate(q, GATEID_Z, 18);    	q = qdd_gate(q, GATEID_Z, 24);    	q = qdd_cgate(q, GATEID_Z, 20, 22);	q = qdd_gate(q, GATEID_S, 3);    
    q = qdd_cgate(q, GATEID_X, 1, 3);	q = qdd_cgate(q, GATEID_Z, 2, 6);	q = qdd_gate(q, GATEID_H, 22);    	q = qdd_gate(q, GATEID_H, 3);    
    q = qdd_gate(q, GATEID_Z, 23);    	q = qdd_gate(q, GATEID_S, 6);    	q = qdd_gate(q, GATEID_S, 24);    	q = qdd_cgate(q, GATEID_Z, 2, 17);
    q = qdd_cgate(q, GATEID_Z, 1, 4);	q = qdd_gate(q, GATEID_H, 23);    	q = qdd_cgate(q, GATEID_Z, 0, 10);	q = qdd_cgate(q, GATEID_Z, 1, 24);
    q = qdd_gate(q, GATEID_X, 10);    	q = qdd_cgate(q, GATEID_Z, 15, 21);	q = qdd_cgate(q, GATEID_X, 2, 8);	q = qdd_gate(q, GATEID_H, 18);    
    q = qdd_cgate(q, GATEID_X, 1, 20);	q = qdd_cgate(q, GATEID_X, 14, 23);	q = qdd_gate(q, GATEID_Z, 8);    	q = qdd_gate(q, GATEID_X, 10);    
    q = qdd_gate(q, GATEID_H, 14);    	q = qdd_cgate(q, GATEID_X, 19, 20);	q = qdd_cgate(q, GATEID_Z, 22, 23);	q = qdd_cgate(q, GATEID_Z, 0, 8);
    q = qdd_cgate(q, GATEID_X, 9, 14);	q = qdd_cgate(q, GATEID_X, 10, 14);	q = qdd_cgate(q, GATEID_Z, 1, 21);	q = qdd_cgate(q, GATEID_X, 5, 13);
    q = qdd_cgate(q, GATEID_Z, 6, 9);	q = qdd_cgate(q, GATEID_Z, 16, 22);	q = qdd_cgate(q, GATEID_X, 15, 17);	q = qdd_cgate(q, GATEID_X, 13, 14);
    q = qdd_cgate(q, GATEID_X, 4, 10);	q = qdd_gate(q, GATEID_H, 19);    	q = qdd_cgate(q, GATEID_Z, 6, 24);	q = qdd_cgate(q, GATEID_Z, 17, 23);
    q = qdd_cgate(q, GATEID_X, 6, 13);	q = qdd_cgate(q, GATEID_X, 9, 24);	q = qdd_cgate(q, GATEID_X, 6, 12);	q = qdd_cgate(q, GATEID_X, 7, 15);
    q = qdd_cgate(q, GATEID_Z, 17, 23);	q = qdd_cgate(q, GATEID_X, 14, 21);	q = qdd_gate(q, GATEID_H, 3);    	q = qdd_gate(q, GATEID_S, 0);    
    q = qdd_cgate(q, GATEID_Z, 10, 20);	q = qdd_cgate(q, GATEID_X, 5, 9);	q = qdd_gate(q, GATEID_H, 19);    	q = qdd_gate(q, GATEID_H, 1);    
    q = qdd_cgate(q, GATEID_Z, 5, 22);	q = qdd_cgate(q, GATEID_Z, 12, 21);	q = qdd_gate(q, GATEID_S, 5);    	q = qdd_gate(q, GATEID_H, 22);    
    q = qdd_cgate(q, GATEID_Z, 0, 21);	q = qdd_gate(q, GATEID_Z, 17);    	q = qdd_cgate(q, GATEID_Z, 16, 23);	q = qdd_gate(q, GATEID_H, 15);    
    q = qdd_cgate(q, GATEID_X, 5, 18);	q = qdd_gate(q, GATEID_X, 9);    	q = qdd_gate(q, GATEID_S, 1);    	q = qdd_cgate(q, GATEID_Z, 15, 21);
    q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_Z, 15);    	q = qdd_gate(q, GATEID_X, 20);    	q = qdd_cgate(q, GATEID_Z, 7, 9);
    q = qdd_gate(q, GATEID_Z, 18);    	q = qdd_gate(q, GATEID_X, 18);    	q = qdd_cgate(q, GATEID_Z, 15, 19);	q = qdd_cgate(q, GATEID_Z, 17, 24);
    q = qdd_gate(q, GATEID_S, 17);    	q = qdd_cgate(q, GATEID_X, 14, 20);	q = qdd_gate(q, GATEID_S, 7);    	q = qdd_cgate(q, GATEID_Z, 10, 19);
    q = qdd_cgate(q, GATEID_Z, 0, 11);	q = qdd_gate(q, GATEID_Z, 15);    	q = qdd_gate(q, GATEID_X, 0);    	q = qdd_gate(q, GATEID_H, 7);    
    q = qdd_gate(q, GATEID_H, 7);    	q = qdd_cgate(q, GATEID_Z, 10, 11);	q = qdd_cgate(q, GATEID_X, 1, 5);	q = qdd_cgate(q, GATEID_Z, 7, 15);
    q = qdd_cgate(q, GATEID_X, 1, 17);	q = qdd_cgate(q, GATEID_Z, 2, 16);	q = qdd_gate(q, GATEID_S, 4);    	q = qdd_cgate(q, GATEID_X, 9, 18);
    q = qdd_gate(q, GATEID_H, 17);    	q = qdd_cgate(q, GATEID_Z, 1, 2);	q = qdd_cgate(q, GATEID_X, 7, 10);	q = qdd_cgate(q, GATEID_Z, 2, 20);
    q = qdd_cgate(q, GATEID_Z, 10, 18);	q = qdd_cgate(q, GATEID_X, 3, 15);	q = qdd_gate(q, GATEID_S, 1);    	q = qdd_cgate(q, GATEID_X, 4, 5);
    q = qdd_gate(q, GATEID_H, 18);    	q = qdd_cgate(q, GATEID_X, 4, 13);	q = qdd_gate(q, GATEID_H, 22);    	q = qdd_gate(q, GATEID_Z, 8);    
    q = qdd_cgate(q, GATEID_X, 12, 23);	q = qdd_gate(q, GATEID_H, 22);    	q = qdd_cgate(q, GATEID_Z, 15, 21);	q = qdd_cgate(q, GATEID_X, 5, 18);
    q = qdd_cgate(q, GATEID_Z, 4, 18);	q = qdd_gate(q, GATEID_H, 22);    	q = qdd_cgate(q, GATEID_X, 11, 15);	q = qdd_gate(q, GATEID_H, 10);    
    q = qdd_gate(q, GATEID_H, 7);    	q = qdd_cgate(q, GATEID_Z, 9, 21);	q = qdd_cgate(q, GATEID_Z, 21, 24);	q = qdd_gate(q, GATEID_X, 14);    
    q = qdd_gate(q, GATEID_X, 8);    	q = qdd_cgate(q, GATEID_Z, 8, 14);	q = qdd_cgate(q, GATEID_Z, 10, 17);	q = qdd_gate(q, GATEID_X, 24);    
    q = qdd_cgate(q, GATEID_X, 4, 24);	q = qdd_cgate(q, GATEID_X, 9, 15);	q = qdd_gate(q, GATEID_X, 9);    	q = qdd_gate(q, GATEID_S, 20);    
    q = qdd_cgate(q, GATEID_Z, 10, 23);	q = qdd_gate(q, GATEID_H, 16);    	q = qdd_cgate(q, GATEID_Z, 6, 9);	q = qdd_gate(q, GATEID_Z, 23);    
    q = qdd_gate(q, GATEID_H, 24);    	q = qdd_cgate(q, GATEID_Z, 12, 21);	q = qdd_gate(q, GATEID_X, 4);    	q = qdd_cgate(q, GATEID_X, 0, 11);
    q = qdd_gate(q, GATEID_H, 20);    	q = qdd_gate(q, GATEID_H, 18);    	q = qdd_gate(q, GATEID_S, 19);    	q = qdd_cgate(q, GATEID_Z, 6, 13);
    q = qdd_cgate(q, GATEID_X, 21, 23);	q = qdd_cgate(q, GATEID_Z, 13, 23);	q = qdd_gate(q, GATEID_H, 6);    	q = qdd_gate(q, GATEID_Z, 14);    
    q = qdd_cgate(q, GATEID_Z, 10, 23);	q = qdd_gate(q, GATEID_S, 15);    	q = qdd_gate(q, GATEID_H, 9);    	q = qdd_gate(q, GATEID_S, 20);    
    q = qdd_cgate(q, GATEID_X, 2, 19);	q = qdd_gate(q, GATEID_H, 16);    	q = qdd_cgate(q, GATEID_Z, 6, 20);	q = qdd_cgate(q, GATEID_Z, 9, 22);
    q = qdd_gate(q, GATEID_H, 8);    	q = qdd_cgate(q, GATEID_Z, 6, 7);	q = qdd_cgate(q, GATEID_X, 3, 24);	q = qdd_cgate(q, GATEID_X, 6, 22);
    q = qdd_cgate(q, GATEID_X, 6, 7);	q = qdd_cgate(q, GATEID_X, 4, 17);	q = qdd_gate(q, GATEID_S, 10);    	q = qdd_gate(q, GATEID_S, 18);    
    q = qdd_gate(q, GATEID_S, 3);    	q = qdd_gate(q, GATEID_X, 4);    	q = qdd_cgate(q, GATEID_X, 14, 21);	q = qdd_cgate(q, GATEID_X, 19, 20);
    q = qdd_gate(q, GATEID_Z, 21);    	q = qdd_gate(q, GATEID_H, 10);    	q = qdd_cgate(q, GATEID_Z, 18, 24);	q = qdd_cgate(q, GATEID_Z, 4, 13);
    q = qdd_gate(q, GATEID_S, 23);    	q = qdd_gate(q, GATEID_Z, 10);    	q = qdd_cgate(q, GATEID_Z, 2, 7);	q = qdd_cgate(q, GATEID_X, 0, 21);
    q = qdd_cgate(q, GATEID_Z, 2, 23);	q = qdd_cgate(q, GATEID_X, 10, 18);	q = qdd_cgate(q, GATEID_Z, 10, 19);	q = qdd_cgate(q, GATEID_Z, 1, 23);
    q = qdd_cgate(q, GATEID_Z, 11, 14);	q = qdd_gate(q, GATEID_X, 6);    	q = qdd_gate(q, GATEID_H, 4);    	q = qdd_cgate(q, GATEID_Z, 6, 24);
    q = qdd_cgate(q, GATEID_X, 13, 18);	q = qdd_gate(q, GATEID_S, 11);    	q = qdd_cgate(q, GATEID_Z, 14, 17);	q = qdd_gate(q, GATEID_X, 23);    
    q = qdd_cgate(q, GATEID_Z, 7, 23);	q = qdd_gate(q, GATEID_Z, 2);    	q = qdd_cgate(q, GATEID_Z, 2, 11);	q = qdd_cgate(q, GATEID_X, 0, 22);
    q = qdd_cgate(q, GATEID_Z, 1, 6);	q = qdd_cgate(q, GATEID_X, 12, 13);	q = qdd_cgate(q, GATEID_Z, 6, 13);	q = qdd_gate(q, GATEID_S, 6);    
    q = qdd_gate(q, GATEID_H, 4);    	q = qdd_gate(q, GATEID_Z, 2);    	q = qdd_cgate(q, GATEID_X, 2, 5);	q = qdd_gate(q, GATEID_H, 17);    
    q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_X, 23);    	q = qdd_gate(q, GATEID_Z, 18);    	q = qdd_gate(q, GATEID_S, 15);    
    q = qdd_gate(q, GATEID_S, 13);    	q = qdd_cgate(q, GATEID_Z, 14, 22);	q = qdd_gate(q, GATEID_H, 4);    	q = qdd_cgate(q, GATEID_X, 8, 19);
    q = qdd_cgate(q, GATEID_X, 3, 13);	q = qdd_gate(q, GATEID_H, 10);    	q = qdd_cgate(q, GATEID_Z, 12, 24);	q = qdd_gate(q, GATEID_X, 15);    
    q = qdd_cgate(q, GATEID_Z, 1, 19);	q = qdd_gate(q, GATEID_H, 24);    	q = qdd_gate(q, GATEID_S, 15);    	q = qdd_gate(q, GATEID_S, 7);    
    q = qdd_gate(q, GATEID_S, 3);    	q = qdd_gate(q, GATEID_X, 22);    	q = qdd_gate(q, GATEID_H, 7);    	q = qdd_cgate(q, GATEID_X, 3, 14);
    q = qdd_cgate(q, GATEID_Z, 4, 13);	q = qdd_cgate(q, GATEID_X, 12, 23);	q = qdd_gate(q, GATEID_X, 17);    	q = qdd_gate(q, GATEID_S, 24);    
    q = qdd_cgate(q, GATEID_X, 21, 23);	q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_X, 0);    	q = qdd_gate(q, GATEID_H, 6);    
    q = qdd_gate(q, GATEID_H, 22);    	q = qdd_cgate(q, GATEID_Z, 14, 16);	q = qdd_cgate(q, GATEID_X, 2, 17);	q = qdd_cgate(q, GATEID_Z, 1, 17);
    q = qdd_cgate(q, GATEID_Z, 1, 16);	q = qdd_gate(q, GATEID_Z, 13);    	q = qdd_cgate(q, GATEID_Z, 6, 18);	q = qdd_cgate(q, GATEID_Z, 3, 12);
    q = qdd_cgate(q, GATEID_X, 5, 13);	q = qdd_cgate(q, GATEID_Z, 13, 17);	q = qdd_gate(q, GATEID_X, 6);    	q = qdd_gate(q, GATEID_H, 21);    
    q = qdd_cgate(q, GATEID_X, 1, 20);	q = qdd_gate(q, GATEID_Z, 6);    	q = qdd_gate(q, GATEID_H, 5);    	q = qdd_cgate(q, GATEID_Z, 7, 21);
    q = qdd_cgate(q, GATEID_Z, 7, 10);	q = qdd_cgate(q, GATEID_X, 7, 12);	q = qdd_cgate(q, GATEID_X, 4, 19);	q = qdd_cgate(q, GATEID_Z, 1, 14);
    q = qdd_gate(q, GATEID_H, 24);    	q = qdd_gate(q, GATEID_X, 10);    	q = qdd_cgate(q, GATEID_X, 2, 16);	q = qdd_cgate(q, GATEID_X, 2, 7);
    q = qdd_gate(q, GATEID_S, 6);    	q = qdd_gate(q, GATEID_X, 8);    	q = qdd_cgate(q, GATEID_Z, 2, 16);	q = qdd_gate(q, GATEID_Z, 14);    
    q = qdd_cgate(q, GATEID_X, 4, 8);	q = qdd_gate(q, GATEID_S, 2);    	q = qdd_cgate(q, GATEID_Z, 11, 20);	q = qdd_gate(q, GATEID_Z, 16);    
    q = qdd_gate(q, GATEID_X, 22);    	q = qdd_cgate(q, GATEID_X, 8, 17);	q = qdd_cgate(q, GATEID_X, 0, 3);	q = qdd_gate(q, GATEID_H, 22);    
    q = qdd_gate(q, GATEID_X, 4);    	q = qdd_gate(q, GATEID_H, 1);    	q = qdd_cgate(q, GATEID_X, 11, 22);	q = qdd_cgate(q, GATEID_Z, 10, 18);
    q = qdd_gate(q, GATEID_S, 20);    	q = qdd_cgate(q, GATEID_X, 5, 9);	q = qdd_cgate(q, GATEID_X, 8, 11);	q = qdd_cgate(q, GATEID_X, 1, 16);
    q = qdd_cgate(q, GATEID_Z, 17, 18);	q = qdd_cgate(q, GATEID_X, 3, 8);	q = qdd_gate(q, GATEID_S, 18);    	q = qdd_cgate(q, GATEID_X, 4, 19);
    q = qdd_gate(q, GATEID_Z, 19);    	q = qdd_gate(q, GATEID_H, 0);    	q = qdd_cgate(q, GATEID_X, 3, 24);	q = qdd_cgate(q, GATEID_Z, 8, 16);
    q = qdd_cgate(q, GATEID_X, 0, 24);	q = qdd_cgate(q, GATEID_X, 13, 21);	q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_H, 13);    
    q = qdd_gate(q, GATEID_Z, 20);    	q = qdd_cgate(q, GATEID_Z, 0, 5);	q = qdd_gate(q, GATEID_H, 4);    	q = qdd_gate(q, GATEID_Z, 15);    
    q = qdd_cgate(q, GATEID_Z, 10, 24);	q = qdd_gate(q, GATEID_H, 22);    	q = qdd_cgate(q, GATEID_Z, 21, 23);	q = qdd_gate(q, GATEID_Z, 0);    
    q = qdd_gate(q, GATEID_Z, 14);    	q = qdd_cgate(q, GATEID_X, 8, 17);	q = qdd_gate(q, GATEID_H, 5);    	q = qdd_gate(q, GATEID_X, 11);    
    q = qdd_cgate(q, GATEID_X, 15, 20);	q = qdd_gate(q, GATEID_X, 22);    	q = qdd_gate(q, GATEID_S, 17);    	q = qdd_cgate(q, GATEID_X, 5, 21);
    q = qdd_gate(q, GATEID_S, 23);    	q = qdd_gate(q, GATEID_Z, 13);    	q = qdd_gate(q, GATEID_X, 13);    	q = qdd_cgate(q, GATEID_Z, 7, 12);
    q = qdd_gate(q, GATEID_H, 2);    	q = qdd_gate(q, GATEID_Z, 17);    	q = qdd_cgate(q, GATEID_X, 4, 24);	q = qdd_cgate(q, GATEID_Z, 19, 22);
    q = qdd_gate(q, GATEID_H, 1);    	q = qdd_cgate(q, GATEID_Z, 0, 16);	q = qdd_cgate(q, GATEID_Z, 2, 12);	q = qdd_cgate(q, GATEID_Z, 2, 11);
    q = qdd_cgate(q, GATEID_X, 17, 24);	q = qdd_cgate(q, GATEID_Z, 13, 20);	q = qdd_cgate(q, GATEID_Z, 1, 22);	q = qdd_gate(q, GATEID_S, 24);    
    q = qdd_gate(q, GATEID_Z, 24);    	q = qdd_gate(q, GATEID_H, 23);    	q = qdd_cgate(q, GATEID_X, 2, 16);	q = qdd_cgate(q, GATEID_Z, 2, 4);
    q = qdd_gate(q, GATEID_X, 17);    	q = qdd_cgate(q, GATEID_Z, 2, 16);	q = qdd_gate(q, GATEID_S, 11);    	q = qdd_gate(q, GATEID_H, 0);    
    q = qdd_cgate(q, GATEID_Z, 8, 20);	q = qdd_gate(q, GATEID_H, 6);    	q = qdd_cgate(q, GATEID_Z, 3, 17);	q = qdd_gate(q, GATEID_Z, 13);    
    q = qdd_gate(q, GATEID_S, 12);    	q = qdd_gate(q, GATEID_H, 23);    	q = qdd_gate(q, GATEID_S, 1);    	q = qdd_cgate(q, GATEID_Z, 3, 17);
    q = qdd_gate(q, GATEID_Z, 11);    	q = qdd_gate(q, GATEID_H, 24);    	q = qdd_cgate(q, GATEID_Z, 0, 18);	q = qdd_gate(q, GATEID_S, 23);    
    q = qdd_gate(q, GATEID_H, 4);    	q = qdd_gate(q, GATEID_H, 24);    	q = qdd_cgate(q, GATEID_X, 0, 17);	q = qdd_cgate(q, GATEID_X, 16, 23);
    q = qdd_cgate(q, GATEID_Z, 11, 16);	q = qdd_cgate(q, GATEID_X, 7, 9);	q = qdd_cgate(q, GATEID_X, 10, 23);	q = qdd_gate(q, GATEID_H, 11);    
    q = qdd_cgate(q, GATEID_X, 7, 12);	q = qdd_gate(q, GATEID_X, 4);    	q = qdd_cgate(q, GATEID_Z, 5, 14);	q = qdd_cgate(q, GATEID_Z, 6, 24);
    q = qdd_gate(q, GATEID_H, 13);    	q = qdd_gate(q, GATEID_H, 17);    	q = qdd_cgate(q, GATEID_Z, 18, 22);	q = qdd_cgate(q, GATEID_Z, 18, 23);
    q = qdd_cgate(q, GATEID_Z, 18, 21);	q = qdd_cgate(q, GATEID_Z, 2, 8);	q = qdd_gate(q, GATEID_Z, 7);    	q = qdd_cgate(q, GATEID_X, 7, 22);
    q = qdd_cgate(q, GATEID_X, 9, 14);	q = qdd_gate(q, GATEID_H, 17);    	q = qdd_cgate(q, GATEID_Z, 6, 22);	q = qdd_gate(q, GATEID_Z, 20);    
    q = qdd_cgate(q, GATEID_Z, 15, 16);	q = qdd_cgate(q, GATEID_Z, 11, 17);	q = qdd_cgate(q, GATEID_Z, 13, 15);	q = qdd_gate(q, GATEID_H, 3);    
    q = qdd_cgate(q, GATEID_Z, 5, 18);	q = qdd_cgate(q, GATEID_X, 16, 24);	q = qdd_cgate(q, GATEID_X, 0, 20);	q = qdd_gate(q, GATEID_H, 12);    
    q = qdd_gate(q, GATEID_X, 8);    	q = qdd_cgate(q, GATEID_Z, 1, 4);	q = qdd_cgate(q, GATEID_X, 6, 24);	q = qdd_cgate(q, GATEID_X, 14, 18);
    q = qdd_cgate(q, GATEID_X, 6, 10);	q = qdd_cgate(q, GATEID_Z, 7, 15);	q = qdd_gate(q, GATEID_H, 18);    	q = qdd_gate(q, GATEID_X, 23);    
    q = qdd_gate(q, GATEID_S, 10);    	q = qdd_cgate(q, GATEID_X, 18, 24);	q = qdd_gate(q, GATEID_S, 16);    	q = qdd_cgate(q, GATEID_Z, 0, 13);
    q = qdd_cgate(q, GATEID_Z, 5, 6);	q = qdd_cgate(q, GATEID_Z, 5, 23);	q = qdd_gate(q, GATEID_S, 8);    	q = qdd_gate(q, GATEID_S, 2);    
    q = qdd_gate(q, GATEID_X, 7);    	q = qdd_gate(q, GATEID_Z, 7);    	q = qdd_cgate(q, GATEID_X, 3, 7);	q = qdd_gate(q, GATEID_H, 23);    
    q = qdd_gate(q, GATEID_Z, 16);    	q = qdd_cgate(q, GATEID_Z, 0, 11);	q = qdd_cgate(q, GATEID_Z, 3, 10);	q = qdd_gate(q, GATEID_S, 16);    
    q = qdd_gate(q, GATEID_X, 18);    	q = qdd_cgate(q, GATEID_X, 16, 23);	q = qdd_gate(q, GATEID_S, 19);    	q = qdd_cgate(q, GATEID_Z, 7, 11);
    q = qdd_gate(q, GATEID_H, 18);    	q = qdd_cgate(q, GATEID_X, 11, 22);	q = qdd_gate(q, GATEID_H, 15);    	q = qdd_gate(q, GATEID_H, 0);    
    q = qdd_cgate(q, GATEID_Z, 20, 24);	q = qdd_cgate(q, GATEID_X, 9, 14);	q = qdd_cgate(q, GATEID_Z, 8, 18);	q = qdd_cgate(q, GATEID_X, 2, 14);
    q = qdd_cgate(q, GATEID_Z, 15, 20);	q = qdd_cgate(q, GATEID_X, 19, 24);	q = qdd_cgate(q, GATEID_X, 2, 20);	q = qdd_gate(q, GATEID_H, 15);    
    q = qdd_cgate(q, GATEID_Z, 2, 13);	q = qdd_gate(q, GATEID_H, 0);    	q = qdd_gate(q, GATEID_S, 13);    	q = qdd_gate(q, GATEID_H, 24);    
    
    t_end = wctime();
    runtime = (t_end - t_start);

    node_count = qdd_countnodes(q);
    printf("%ld nodes, %lf sec\n", node_count, runtime);

    // Cleanup
    sylvan_quit();
    lace_exit();
    
    return 0;
}


int bench_supremacy_5_1(uint32_t depth, uint32_t workers)
{
    // 5x1 "grid" from [Characterizing Quantum Supremacy in Near-Term Devices]
    printf("bench sup5_1, depth %3d, %2d worker(s), ", depth, workers);
    fflush(stdout);

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
    sylvan_init_qdd(1LL<<23, -1);

    srand(42);
    QDD res = supremacy_5_1_circuit(depth);
    node_count = qdd_countnodes(res);

    t_end = wctime();
    runtime = (t_end - t_start);

    printf("%4ld nodes, %lf sec\n", node_count, runtime);

    // Cleanup
    sylvan_quit();
    lace_exit();

    return 0;
}

double bench_supremacy_5_4_once(uint32_t depth, uint32_t workers, uint64_t rseed, char *fpath, uint64_t *nodes_peak, uint64_t *n_gates)
{
    if (VERBOSE) {
        printf("bench sup5_4, depth %3d, %2d worker(s), ", depth, workers);
        fflush(stdout);
    }

    FILE *logfile = NULL;
    if (fpath != NULL)
        logfile = fopen(fpath, "w");

    uint64_t node_count;
    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    sylvan_init_qdd(ctable_size, ctable_tolerance);
    qdd_set_gc_ctable_thres(ctable_gc_thres);

    qdd_stats_start(logfile);

    srand(rseed);
    QDD res = supremacy_5_4_circuit(depth);

    if (nodes_peak != NULL) *nodes_peak = qdd_stats_get_nodes_peak();
    if (n_gates    != NULL) *n_gates = qdd_stats_get_logcounter();
    if (logfile    != NULL) qdd_stats_finish();

    t_end = wctime();
    runtime = (t_end - t_start);

    if (VERBOSE) {
        node_count = qdd_countnodes(res); 
        printf("%4ld nodes, %lf sec\n", node_count, runtime);
    }

    // Cleanup
    sylvan_quit();
    lace_exit();

    return runtime;
}


double bench_grover_once(int num_qubits, bool flag[], int workers, char *fpath, uint64_t *nodes_peak, uint64_t *n_gates)
{
    if (VERBOSE) {
        printf("bench grover, %d qubits, %2d worker(s), ", num_qubits, workers); 
        printf("flag = [");
        for (int i = 0; i < num_qubits; i++)
            printf("%d",flag[i]);
        printf("], ");
        fflush(stdout);
    }

    FILE *logfile = NULL;
    if (fpath != NULL)
        logfile = fopen(fpath, "w");

    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    sylvan_init_qdd(ctable_size, ctable_tolerance);
    qdd_set_gc_ctable_thres(ctable_gc_thres);

    QDD grov;
    uint64_t node_count_end;
    
    qdd_stats_start(logfile);

    grov = qdd_grover(num_qubits, flag);

    if (nodes_peak != NULL) *nodes_peak = qdd_stats_get_nodes_peak();
    if (n_gates    != NULL) *n_gates = qdd_stats_get_logcounter();
    if (logfile    != NULL) qdd_stats_finish();

    t_end = wctime();
    runtime = (t_end - t_start);

    if (VERBOSE) {
        node_count_end = qdd_countnodes(grov);
        double prob = comp_to_prob(comp_value(qdd_get_amplitude(grov, flag)));
        uint64_t np = (nodes_peak == NULL) ? 0 : *nodes_peak;
        printf("%ld nodes end (%ld peak), Pr(flag)=%.3lf, %lf sec\n", node_count_end, np, prob, runtime);
    }

    if (logfile != NULL)
        fclose(logfile);

    // Cleanup
    sylvan_quit();
    lace_exit();
    return runtime;
}

double bench_shor_once(uint64_t N, uint64_t a, int workers, int rseed, bool *success, char *fpath, uint64_t *nodes_peak, uint64_t *n_gates)
{
    if (VERBOSE) {
        uint32_t num_qubits = (int)ceil(log2(N))*2 + 3;
        printf("bench shor, factor %ld (%d qubits), %2d worker(s), ", N, num_qubits, workers); 
        fflush(stdout);
    }

    FILE *logfile = NULL;
    if (fpath != NULL)
        logfile = fopen(fpath, "w");

    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    sylvan_init_qdd(ctable_size, ctable_tolerance);

    qdd_stats_start(logfile);

    srand(rseed);
    uint64_t fac = run_shor(N, a, false);

    t_end = wctime();
    runtime = (t_end - t_start);

    if (nodes_peak != NULL) *nodes_peak = qdd_stats_get_nodes_peak();
    if (n_gates    != NULL) *n_gates = qdd_stats_get_logcounter();
    if (logfile    != NULL) {
        qdd_stats_finish();
        fclose(logfile);
    }
    *success = (fac == 0) ? 0 : 1;

    if (VERBOSE) printf("found factor %ld, %lf sec\n", fac, runtime);

    // Cleanup
    sylvan_quit();
    lace_exit();
    return runtime;
}

int bench_supremacy()
{
    VERBOSE = true;

    // output dir
    mkdir("benchmark_data/supremacy/", 0700);
    char output_dir[256];
    sprintf(output_dir, "benchmark_data/supremacy/%ld/", time(NULL));
    mkdir(output_dir, 0700);
    char history_dir[256];
    strcpy(history_dir, output_dir);
    strcat(history_dir, "run_histories/");
    mkdir(history_dir, 0700);
    // output file for runtime data
    char overview_fname[256];
    strcpy(overview_fname, output_dir);
    strcat(overview_fname, "summary.csv");
    FILE *overview_file = fopen(overview_fname, "w");
    fprintf(overview_file, "qubits, depth, rseed, peak_nodes, workers, "
                           "gates, runtime, avg_gate_time, "
                           "plus_cacheput, plus_cached, "
                           "gate_cacheput, gate_cached, "
                           "cgate_cacheput, cgate_cached\n");
    // output file for sylvan parameters
    char param_fname[256];
    strcpy(param_fname, output_dir);
    strcat(param_fname, "parameters.json");
    FILE *param_file = fopen(param_fname, "w");

    // sylvan / qdd params
    min_tablesize = max_tablesize = 1LL<<30;
    min_cachesize = max_cachesize = 1LL<<16;
    ctable_size   = 1LL<<23;
    ctable_gc_thres = 0.25;
    ctable_tolerance = 1e-14;
    write_parameters(param_file);

    // params
    int nqubits = 20; // always 20 for 5x4 grid
    int depths[] = {15};//{15,16,17,18,19,20};
    int ndepths = 1;

    // different number of workers to test
    int n_workers[] = {1, 2, 4};
    int nn_workers  = 3;

    // re-runs for different depths
    int re_runs = 2;
    uint64_t rseeds[] = {66, 123};

    // runtimes are written to single file
    double runtime, avg_gate_time;
    uint64_t nodes_peak, n_gates;
    uint64_t plus_cacheput, gate_cacheput, cgate_cacheput;
    uint64_t plus_cached, gate_cached, cgate_cached;

    for (int i = 0; i < ndepths; i++) {
        for (int r = 0; r < re_runs; r++) {
            uint64_t rseed = rseeds[r];
            for (int w = 0; w < nn_workers; w++) {

                // output file for history of this run
                char history_path[256];
                char history_fname[256];
                sprintf(history_fname, "sup5x4_d%d_w%d_rseed%ld.csv", depths[i], n_workers[w], rseed);
                strcpy(history_path, history_dir);
                strcat(history_path, history_fname);

                // bench twice, once with logging and once for timing
                runtime = bench_supremacy_5_4_once(depths[i], n_workers[w], rseed, NULL, NULL, NULL);
                bench_supremacy_5_4_once(depths[i], n_workers[w], rseed, history_path, &nodes_peak, &n_gates);

                // add summary of this run to overview file
                avg_gate_time = runtime / (double) n_gates;
                #if SYLVAN_STATS
                plus_cacheput  = sylvan_stats.counters[QDD_PLUS_CACHEDPUT];
                gate_cacheput  = sylvan_stats.counters[QDD_GATE_CACHEDPUT];
                cgate_cacheput = sylvan_stats.counters[QDD_CGATE_CACHEDPUT];
                plus_cached    = sylvan_stats.counters[QDD_PLUS_CACHED];
                gate_cached    = sylvan_stats.counters[QDD_GATE_CACHED];
                cgate_cached   = sylvan_stats.counters[QDD_CGATE_CACHED];
                #else
                plus_cached = gate_cached = cgate_cached = 0;
                plus_cacheput = gate_cacheput = cgate_cacheput = 0;
                #endif
                fprintf(overview_file, "%d, %d, %ld, %ld, %d, %ld, %lf, %.3e, %ld, %ld, %ld, %ld, %ld, %ld\n",
                                        nqubits, depths[i], rseed, nodes_peak, n_workers[w],
                                        n_gates, runtime, avg_gate_time,
                                        plus_cacheput, plus_cached,
                                        gate_cacheput, gate_cached,
                                        cgate_cacheput, cgate_cached);
            }
        }
    }
    

    return 0;
}

int bench_grover()
{
    VERBOSE = true;
    
    // output dir
    mkdir("benchmark_data/grover/", 0700);
    char output_dir[256];
    sprintf(output_dir, "benchmark_data/grover/%ld/", time(NULL));
    mkdir(output_dir, 0700);
    char history_dir[256];
    strcpy(history_dir, output_dir);
    strcat(history_dir, "run_histories/");
    mkdir(history_dir, 0700);
    // output file for runtime data
    char overview_fname[256];
    strcpy(overview_fname, output_dir);
    strcat(overview_fname, "summary.csv");
    FILE *overview_file = fopen(overview_fname, "w");
    fprintf(overview_file, "qubits, peak_nodes, workers, "
                           "gates, runtime, avg_gate_time, "
                           "plus_cacheput, plus_cached, "
                           "gate_cacheput, gate_cached, "
                           "cgate_cacheput, cgate_cached, "
                           "flag\n");
    // output file for sylvan parameters
    char param_fname[256];
    strcpy(param_fname, output_dir);
    strcat(param_fname, "parameters.json");
    FILE *param_file = fopen(param_fname, "w");

    // sylvan / qdd params
    min_tablesize = max_tablesize = 1LL<<25;
    min_cachesize = max_cachesize = 1LL<<16;
    ctable_size   = 1LL<<18;
    ctable_tolerance = 1e-14;
    write_parameters(param_file);

    // different number of qubits to test
    int n_qubits[] = {15, 20};
    int nn_qubits  = 2;
    
    // different number of workers to test
    int n_workers[] = {1, 2, 4};
    int nn_workers  = 3;

    // different number of random flags to test
    int n_flags = 3;
    bool *flag;
    int f_int;

    // runtimes are written to single file
    double runtime, avg_gate_time;
    uint64_t nodes_peak, n_gates;
    uint64_t plus_cacheput, gate_cacheput, cgate_cacheput;
    uint64_t plus_cached, gate_cached, cgate_cached;

    // run benchmarks
    srand(42);
    for (int q = 0; q < nn_qubits; q++) {

        for (int f = 0; f < n_flags; f++) {
            flag  = qdd_grover_random_flag(n_qubits[q]);
            f_int = bitarray_to_int(flag, n_qubits[q], true);

            for (int w = 0; w < nn_workers; w++) {

                // output file for history of this run
                char history_path[256];
                char history_fname[256];
                sprintf(history_fname, "grov_hist_n%d_w%d_f%d.csv", n_qubits[q], n_workers[w], f_int);
                strcpy(history_path, history_dir);
                strcat(history_path, history_fname);

                // bench twice, once with logging and once for timing
                runtime = bench_grover_once(n_qubits[q], flag, n_workers[w], NULL, NULL, NULL);
                bench_grover_once(n_qubits[q], flag, n_workers[w], history_path, &nodes_peak, &n_gates);

                // add summary of this run to overview file
                avg_gate_time = runtime / (double) n_gates;
                #if SYLVAN_STATS
                plus_cacheput  = sylvan_stats.counters[QDD_PLUS_CACHEDPUT];
                gate_cacheput  = sylvan_stats.counters[QDD_GATE_CACHEDPUT];
                cgate_cacheput = sylvan_stats.counters[QDD_CGATE_CACHEDPUT];
                plus_cached    = sylvan_stats.counters[QDD_PLUS_CACHED];
                gate_cached    = sylvan_stats.counters[QDD_GATE_CACHED];
                cgate_cached   = sylvan_stats.counters[QDD_CGATE_CACHED];
                #else
                plus_cached = gate_cached = cgate_cached = 0;
                plus_cacheput = gate_cacheput = cgate_cacheput = 0;
                #endif
                fprintf(overview_file, "%d, %ld, %d, %ld, %lf, %.3e, %ld, %ld, %ld, %ld, %ld, %ld, %d\n",
                                        n_qubits[q], nodes_peak, n_workers[w],
                                        n_gates, runtime, avg_gate_time, 
                                        plus_cacheput, plus_cached,
                                        gate_cacheput, gate_cached,
                                        cgate_cacheput, cgate_cached,
                                        f_int);
            }
        }
    }

    fclose(overview_file);

    return 0;
}

int bench_shor()
{
    VERBOSE = true;
    
    // output dir
    mkdir("benchmark_data/shor/", 0700);
    char output_dir[256];
    sprintf(output_dir, "benchmark_data/shor/%ld/", time(NULL));
    mkdir(output_dir, 0700);
    char history_dir[256];
    strcpy(history_dir, output_dir);
    strcat(history_dir, "run_histories/");
    mkdir(history_dir, 0700);
    // output file for runtime data
    char overview_fname[256];
    strcpy(overview_fname, output_dir);
    strcat(overview_fname, "summary.csv");
    FILE *overview_file = fopen(overview_fname, "w");
    fprintf(overview_file, "N, a, qubits, peak_nodes, success, "
                           "workers, gates, runtime, avg_gate_time, "
                           "plus_cacheput, plus_cached, "
                           "gate_cacheput, gate_cached, "
                           "cgate_cacheput, cgate_cached\n");
    // output file for sylvan parameters
    char param_fname[256];
    strcpy(param_fname, output_dir);
    strcat(param_fname, "parameters.json");
    FILE *param_file = fopen(param_fname, "w");

    // sylvan / qdd params
    min_tablesize = max_tablesize = 1LL<<25;
    min_cachesize = max_cachesize = 1LL<<16;
    ctable_size   = 1LL<<18;
    ctable_tolerance = 1e-14;
    write_parameters(param_file);

    // Different sized N to test
    //  3 x  5 =  15 (11 qubits)
    //  5 x  7 =  35 (15 qubits)
    // 11 x 13 = 143 (19 qubits)
    int Ns[] = {15, 35, 143};
    int nN = 3;
    uint64_t N, a;
    uint32_t nqubits;

    // how often to re-run the same (N,a) (TODO: different 'a' for each run?)
    int re_runs = 2;
    
    // different number of workers to test
    int n_workers[] = {1, 2, 4};
    int nn_workers  = 3;

    // runtimes are written to single file
    double runtime, avg_gate_time;
    uint64_t nodes_peak, n_gates;
    uint64_t plus_cacheput, gate_cacheput, cgate_cacheput;
    uint64_t plus_cached, gate_cached, cgate_cached;

    // run benchmarks
    srand(42);
    for (int q = 0; q < nN; q++) {
        N = Ns[q];
        a = shor_generate_a(N);
        nqubits = (int)ceil(log2(N))*2 + 3;
        for (int r = 0; r < re_runs; r++) {
            uint64_t rseed = rand();
            for (int w = 0; w < nn_workers; w++) {
                
                // output file for history of this run
                char history_path[256];
                char history_fname[256];
                sprintf(history_fname, "shor_hist_n%d_w%d_N%ld_a%ld_r%d.csv", nqubits, n_workers[w], N,a,r);
                strcpy(history_path, history_dir);
                strcat(history_path, history_fname);

                // bench twice, once with logging and once for timing
                bool success;
                runtime = bench_shor_once(N, a, n_workers[w], rseed, &success, NULL, NULL, NULL);
                bench_shor_once(N, a, n_workers[w], rseed, &success, history_path, &nodes_peak, &n_gates);

                // add summary of this run to overview file
                avg_gate_time = runtime / (double) n_gates;
                #if SYLVAN_STATS
                plus_cacheput  = sylvan_stats.counters[QDD_PLUS_CACHEDPUT];
                gate_cacheput  = sylvan_stats.counters[QDD_GATE_CACHEDPUT];
                cgate_cacheput = sylvan_stats.counters[QDD_CGATE_CACHEDPUT];
                plus_cached    = sylvan_stats.counters[QDD_PLUS_CACHED];
                gate_cached    = sylvan_stats.counters[QDD_GATE_CACHED];
                cgate_cached   = sylvan_stats.counters[QDD_CGATE_CACHED];
                #else
                plus_cached = gate_cached = cgate_cached = 0;
                plus_cacheput = gate_cacheput = cgate_cacheput = 0;
                #endif
                fprintf(overview_file, "%ld, %ld, %d, %ld, %d, %d, %ld, %lf, %.3e, %ld, %ld, %ld, %ld, %ld, %ld\n",
                                        N, a, nqubits, nodes_peak, success,
                                        n_workers[w], n_gates, runtime, avg_gate_time, 
                                        plus_cacheput, plus_cached,
                                        gate_cacheput, gate_cached,
                                        cgate_cacheput, cgate_cached);
            }
        }
    }

    fclose(overview_file);

    return 0;
}

int main()
{
    #ifdef HAVE_PROFILER
        // TODO: automate file name depending on which circuit / parameteres
        if (profile_name != NULL) {
            printf("writing profile to %s\n", profile_name);
            ProfilerStart(profile_name);
        }
    #endif

    mkdir("benchmark_data", 0700);
  
    //bench_grover();
    //bench_shor();
    bench_supremacy();

    #ifdef HAVE_PROFILER
        if (profile_name != NULL) ProfilerStop();
    #endif

    return 0;
}