#include <stdio.h>
#include <sys/time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_int.h"

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
    sylvan_init_qdd(1LL<<20);

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
    free_amplitude_table();
    sylvan_quit();
    lace_exit();
    
    return 0;
}

int bench_grover(int num_qubits, bool flag[], int workers)
{
    printf("bench grover, %d qubits, %2d worker(s), ", num_qubits, workers); 
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
    sylvan_init_qdd(1LL<<20);
    sylvan_gc_disable(); // issue with gc, maybe "marked" flag location MTBBD vs QDD


    QDD grov;
    uint64_t node_count;

    grov = qdd_grover(num_qubits, flag);

    t_end = wctime();
    runtime = (t_end - t_start);

    node_count = qdd_countnodes(grov);
    printf("%ld nodes, %lf sec ", node_count, runtime);
    printf("(flag = ");
    for (int i = 0; i < num_qubits; i++)
        printf("%d",flag[i]);
    printf(")\n");

    // Cleanup
    free_amplitude_table();
    sylvan_quit();
    lace_exit();
    return 0;
}

uint32_t random_sqrtXY()
{
    return ( rand() % 2 ) ? GATEID_sqrtX : GATEID_sqrtY;
}

int bench_supremacy_5_1(uint32_t depth)
{
    // 5x1 "grid" from [Characterizing Quantum Supremacy in Near-Term Devices]

    LACE_ME;
    BDDVAR n_qubits = 5;
    srand(42);

    // Start with |00000> state
    QDD qdd = qdd_create_all_zero_state(n_qubits);

    // H on all qubits
    for (BDDVAR k=0; k < n_qubits; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // First 4 cycles the single qubit gates we apply are all T gates
    if (depth > 0) {
        qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);
        qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);
    }
    if (depth > 1) {
        qdd = qdd_cgate(qdd, GATEID_Z, 1, 2);
        qdd = qdd_gate(qdd, GATEID_T, 0);
        qdd = qdd_gate(qdd, GATEID_T, 3);
        qdd = qdd_gate(qdd, GATEID_T, 4);
    }
    if (depth > 2) {
        qdd = qdd_cgate(qdd, GATEID_Z, 2, 3);
        qdd = qdd_gate(qdd, GATEID_T, 1);
    }
    if (depth > 3) {
        qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);
        qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);
        qdd = qdd_gate(qdd, GATEID_T, 2);
    }

    // Following cycles the single qubit gates are random from {sqrt(X), sqrt(Y)}
    for (uint32_t d = 4; d < depth; d++) {
        switch (d % 3) {
        case 0:
            qdd = qdd_cgate(qdd, GATEID_Z, 0, 1);    // CZ(0,1)
            qdd = qdd_cgate(qdd, GATEID_Z, 3, 4);    // CZ(3,4)
            qdd = qdd_gate(qdd, random_sqrtXY(), 2); // random on q2
            break;
        case 1:
            qdd = qdd_cgate(qdd, GATEID_Z, 1, 2);    // CZ(1,2)
            qdd = qdd_gate(qdd, random_sqrtXY(), 0); // random on q0
            qdd = qdd_gate(qdd, random_sqrtXY(), 3); // random on q3
            qdd = qdd_gate(qdd, random_sqrtXY(), 4); // random on q4
            break;
        case 2:
            qdd = qdd_cgate(qdd, GATEID_Z, 2, 3);    // CZ(2,3)
            qdd = qdd_gate(qdd, random_sqrtXY(), 0); // random on q0
            qdd = qdd_gate(qdd, random_sqrtXY(), 1); // random on q1
            qdd = qdd_gate(qdd, random_sqrtXY(), 4); // random on q4
            break;
        default:
            break;
        }
    }

    return 0;
}

int bench_supremacy_5_4(uint32_t depth)
{
    // 5x4 "grid" from [Characterizing Quantum Supremacy in Near-Term Devices]
    // (maybe a lot of hard-coded stuff but oh well...)

    LACE_ME;
    BDDVAR n_qubits = 20;
    srand(66);

    // Start with |00...0> state
    QDD qdd = qdd_create_all_zero_state(n_qubits);

    // H on all qubits
    for (BDDVAR k=0; k < n_qubits; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // qubits which are not involved in a CZ at cycle (d mod 8) 
    // AND had a CZ applied to them in the previous cycle.
    BDDVAR qubits0[6]  = {1,8,10,14,17,19};
    BDDVAR qubits1[8]  = {2,3,5,6,12,13,15,16};
    BDDVAR qubits2[6]  = {0,1,7,10,17,18};
    BDDVAR qubits3[4]  = {6,8,11,13};
    BDDVAR qubits4[4]  = {5,9,10,12};
    BDDVAR qubits5[8]  = {3,4,6,7,13,14,16,17};
    BDDVAR qubits6[4]  = {1,8,12,19};
    BDDVAR qubits7[10] = {0,2,4,5,7,9,11,13,16,18};

    // First 7 cycles the single qubit gates we apply are all T gates
    // (not sure if we also should do random {sqrt(X), sqrt(Y)} on the
    // remaining qubits, the paper is not super clear)
    if (depth > 0) {
        qdd = qdd_cgate(qdd, GATEID_Z,  2,  3);    // CZ(2,3)
        qdd = qdd_cgate(qdd, GATEID_Z,  5,  6);    // CZ(5,6)
        qdd = qdd_cgate(qdd, GATEID_Z, 12, 13);    // CZ(12,13)
        qdd = qdd_cgate(qdd, GATEID_Z, 15, 16);    // CZ(15,16)
    }
    if (depth > 1) {
        qdd = qdd_cgate(qdd, GATEID_Z,  0,  1);    // CZ(0,1)
        qdd = qdd_cgate(qdd, GATEID_Z,  7,  8);    // CZ(7,8)
        qdd = qdd_cgate(qdd, GATEID_Z, 10, 11);    // CZ(10,11)
        qdd = qdd_cgate(qdd, GATEID_Z, 17, 18);    // CZ(17,18)
        BDDVAR qubits[8] = {2,3,5,6,12,13,15,16};
        for (int i = 0; i < 8; i++)
            qdd = qdd_gate(qdd, GATEID_T, qubits[i]);
    }
    if (depth > 2) {
        qdd = qdd_cgate(qdd, GATEID_Z,  6, 11);    // CZ(6,11)
        qdd = qdd_cgate(qdd, GATEID_Z,  8, 13);    // CZ(8,13)
        BDDVAR qubits[8] = {0,1,7,8,10,11,17,18};
        for (int i = 0; i < 8; i++)
            qdd = qdd_gate(qdd, GATEID_T, qubits[i]);
    }
    if (depth > 3) {
        qdd = qdd_cgate(qdd, GATEID_Z,  5, 10);    // CZ(5,10)
        qdd = qdd_cgate(qdd, GATEID_Z,  7, 12);    // CZ(7,12)
        qdd = qdd_cgate(qdd, GATEID_Z,  9, 14);    // CZ(9,14)
    }
    if (depth > 4) {
        qdd = qdd_cgate(qdd, GATEID_Z,  3,  4);    // CZ(3,4)
        qdd = qdd_cgate(qdd, GATEID_Z,  6,  7);    // CZ(6,7)
        qdd = qdd_cgate(qdd, GATEID_Z, 13, 14);    // CZ(13,14)
        qdd = qdd_cgate(qdd, GATEID_Z, 16, 17);    // CZ(16,17)
        qdd = qdd_gate(qdd, GATEID_T,  9);
        qdd = qdd_gate(qdd, GATEID_T, 14);
    }
    if (depth > 5) {
        qdd = qdd_cgate(qdd, GATEID_Z,  1,  2);    // CZ(1,2)
        qdd = qdd_cgate(qdd, GATEID_Z,  8,  9);    // CZ(8,9)
        qdd = qdd_cgate(qdd, GATEID_Z, 11, 12);    // CZ(11,12)
        qdd = qdd_cgate(qdd, GATEID_Z, 18, 19);    // CZ(18,19)
        qdd = qdd_gate(qdd, GATEID_T, 4);
    }
    if (depth > 6) {
        qdd = qdd_cgate(qdd, GATEID_Z,  0,  5);    // CZ(0,5)
        qdd = qdd_cgate(qdd, GATEID_Z,  2,  7);    // CZ(2,7)
        qdd = qdd_cgate(qdd, GATEID_Z,  4,  9);    // CZ(4,9)
        qdd = qdd_cgate(qdd, GATEID_Z, 11, 16);    // CZ(11,16)
        qdd = qdd_cgate(qdd, GATEID_Z, 13, 18);    // CZ(13,18)
        qdd = qdd_gate(qdd, GATEID_T, 19);
    }
    // Following cycles the single qubit gates are random from {sqrt(X), sqrt(Y)}
    for (uint32_t d = 7; d < depth; d++) {
        switch (d % 8) {
        case 0:
            qdd = qdd_cgate(qdd, GATEID_Z,  2,  3);    // CZ(2,3)
            qdd = qdd_cgate(qdd, GATEID_Z,  5,  6);    // CZ(5,6)
            qdd = qdd_cgate(qdd, GATEID_Z, 12, 13);    // CZ(12,13)
            qdd = qdd_cgate(qdd, GATEID_Z, 15, 16);    // CZ(15,16)
            for (int i = 0; i < 6; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits0[i]);
            break;
        case 1:
            qdd = qdd_cgate(qdd, GATEID_Z,  0,  1);    // CZ(0,1)
            qdd = qdd_cgate(qdd, GATEID_Z,  7,  8);    // CZ(7,8)
            qdd = qdd_cgate(qdd, GATEID_Z, 10, 11);    // CZ(10,11)
            qdd = qdd_cgate(qdd, GATEID_Z, 17, 18);    // CZ(17,18)
            for (int i = 0; i < 8; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits1[i]);
            break;
        case 2:
            qdd = qdd_cgate(qdd, GATEID_Z,  6, 11);    // CZ(6,11)
            qdd = qdd_cgate(qdd, GATEID_Z,  8, 13);    // CZ(8,13)
            for (int i = 0; i < 6; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits2[i]);
            break;
        case 3:
            qdd = qdd_cgate(qdd, GATEID_Z,  5, 10);    // CZ(5,10)
            qdd = qdd_cgate(qdd, GATEID_Z,  7, 12);    // CZ(7,12)
            qdd = qdd_cgate(qdd, GATEID_Z,  9, 14);    // CZ(9,14)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits3[i]);
            break;
        case 4:
            qdd = qdd_cgate(qdd, GATEID_Z,  3,  4);    // CZ(3,4)
            qdd = qdd_cgate(qdd, GATEID_Z,  6,  7);    // CZ(6,7)
            qdd = qdd_cgate(qdd, GATEID_Z, 13, 14);    // CZ(13,14)
            qdd = qdd_cgate(qdd, GATEID_Z, 16, 17);    // CZ(16,17)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits4[i]);
            break;
        case 5:
            qdd = qdd_cgate(qdd, GATEID_Z,  1,  2);    // CZ(1,2)
            qdd = qdd_cgate(qdd, GATEID_Z,  8,  9);    // CZ(8,9)
            qdd = qdd_cgate(qdd, GATEID_Z, 11, 12);    // CZ(11,12)
            qdd = qdd_cgate(qdd, GATEID_Z, 18, 19);    // CZ(18,19)
            for (int i = 0; i < 8; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits5[i]);
            break;
        case 6:
            qdd = qdd_cgate(qdd, GATEID_Z,  0,  5);    // CZ(0,5)
            qdd = qdd_cgate(qdd, GATEID_Z,  2,  7);    // CZ(2,7)
            qdd = qdd_cgate(qdd, GATEID_Z,  4,  9);    // CZ(4,9)
            qdd = qdd_cgate(qdd, GATEID_Z, 11, 16);    // CZ(11,16)
            qdd = qdd_cgate(qdd, GATEID_Z, 13, 18);    // CZ(13,18)
            for (int i = 0; i < 4; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits6[i]);
            break;
        case 7:
            qdd = qdd_cgate(qdd, GATEID_Z,  1,  6);    // CZ(1,6)
            qdd = qdd_cgate(qdd, GATEID_Z,  3,  8);    // CZ(3,8)
            qdd = qdd_cgate(qdd, GATEID_Z, 10, 15);    // CZ(10,15)
            qdd = qdd_cgate(qdd, GATEID_Z, 12, 17);    // CZ(12,17)
            qdd = qdd_cgate(qdd, GATEID_Z, 14, 19);    // CZ(14,19)
            for (int i = 0; i < 10; i++) 
                qdd = qdd_gate(qdd, random_sqrtXY(), qubits7[i]);
            break;
        default:
            break;
        }
    }

    return 0;
}


int main()
{
    //bench_25qubit_circuit(1);
    //bench_25qubit_circuit(8);
    //bench_25qubit_circuit(16);

    int n = 22;
    bool flag[n];
    srand(time(NULL));
    for (int i = 0; i < n; i++) flag[i] = (bool)(rand() % 2);
    bench_grover(n, flag, 1);
    bench_grover(n, flag, 8);
    bench_grover(n, flag, 16);

    return 0;
}
