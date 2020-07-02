#include <stdio.h>
#include <time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_int.h"

bool VERBOSE = true;

int bench_25qubit_circuit()
{
    QDD q;
    uint64_t node_count;
    bool x25[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    time_t t_start, t_end, runtime;

    LACE_ME;

    t_start = time(NULL);

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
    t_end = time(NULL);
    node_count = qdd_countnodes(q);


    runtime = (t_end - t_start);
    if(VERBOSE) printf("qdd 25 qubit circuit:     (%ld nodes, %ld sec)\n", node_count, runtime);
    return 0;
}

int bench_grover(int num_qubits)
{
    QDD grov;
    bool flag[num_qubits];
    time_t t_start, t_end, runtime;
    uint64_t node_count;

    // init random flag
    printf("flag: ");
    srand(time(NULL));
    for (int i = 0; i < num_qubits; i++) {
        flag[i] = (bool)(rand() % 2);
        printf("%d",flag[i]);
    }
    printf("\n");

    t_start = time(NULL);
    grov = qdd_grover(num_qubits, flag);
    t_end = time(NULL);
    node_count = qdd_countnodes(grov);

    runtime = (t_end - t_start);
    if(VERBOSE) printf("grover(%d) circuit:     (%ld, %ld sec)\n", num_qubits, node_count, runtime);
    return 0;
}

int runbench()
{
    //if (bench_25qubit_circuit()) return 1;
    if (bench_grover(9)) return 1;

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
    // we also need init_bdd() because some qdd functions 
    // rely on bdd stuff (like cache)
    sylvan_init_bdd();
    // TODO: make sylvan_init_qdd() function and handle stuff there
    init_amplitude_table(19);

    int res = runbench();

    free_amplitude_table();
    sylvan_quit();
    lace_exit();
    return res;
}
