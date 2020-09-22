#include <stdio.h>
#include <sys/time.h>

#include "sylvan.h"
#include "sylvan_qdd_complex.h"
#include "grover.h"
#include "supremacy.h"

#ifdef HAVE_PROFILER
#include <gperftools/profiler.h>
static char* profile_name = NULL; //"bench_qdd.prof";
#endif

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
    sylvan_init_qdd(1LL<<23);

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

int bench_shor(uint64_t N, uint64_t a, int workers, int rand_seed)
{
    uint32_t num_qubits = (int)ceil(log2(N))*2 + 3;
    printf("bench shor, factor %ld (%d qubits), %2d worker(s), ", N, num_qubits, workers); 
    fflush(stdout);

    double t_start, t_end, runtime;
    t_start = wctime();

    // Init Lace
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Init Sylvan
    sylvan_set_limits(4LL<<30, 1, 6);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<23);
    //sylvan_gc_disable(); // issue with gc, maybe "marked" flag location MTBBD vs QDD

    srand(rand_seed);
    uint64_t fac = run_shor(N, a, false);

    t_end = wctime();
    runtime = (t_end - t_start);

    printf("found factor %ld, %lf sec\n", fac, runtime);

    // Cleanup
    sylvan_quit();
    lace_exit();
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

    //bench_25qubit_circuit(1);
    //bench_25qubit_circuit(2);

    int n = 21;
    bool flag[] = {1,1,1,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,0,1};
    bench_grover(n, flag, 1, NULL);


    bench_supremacy_5_1(100, 1);
    bench_supremacy_5_4(5, 1);

    uint64_t rand_seed = time(NULL);
    bench_shor(77, 0, 1, rand_seed);
    bench_shor(77, 0, 2, rand_seed);
    bench_shor(77, 0, 4, rand_seed);

    #ifdef HAVE_PROFILER
        if (profile_name != NULL) ProfilerStop();
    #endif

    return 0;
}