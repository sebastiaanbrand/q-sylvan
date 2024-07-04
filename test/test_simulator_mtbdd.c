/**
 * Copyright 2024 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <inttypes.h>

#include <sylvan.h>
#include "test_assert.h"
#include <sylvan_int.h>
#include <sylvan_mpc.h>

#include "qsylvan_simulator_mtbdd.h"

int
test_create_all_zero_state_double() // TODO: change into complex_t
{
    BDDVAR n = 2;
    MTBDD a = mtbdd_create_all_zero_state_double(n);

    VecArr_t v_arr[(1 << n)];

    mtbdd_to_vector_array(a, n, COLUMN_WISE_MODE, v_arr);

    print_vector_array(v_arr, n);

    test_assert(v_arr[0]==(double)1.0);

    return 0;
}

int
test_create_all_zero_state_complex()
{
    BDDVAR n = 2;
    MTBDD a = mtbdd_create_all_zero_state_mpc(n);

    mpc_ptr v_arr[(1 << n)];

    mtbdd_to_vector_array_mpc(a, n, COLUMN_WISE_MODE, v_arr);

    print_vector_array_mpc(v_arr, n);

    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    test_assert(mpc_compare( (uint64_t)v_arr[0], (uint64_t)mpc_one));

    mpc_clear(mpc_one);

    return 0;
}

int
test_create_basis_state_double() // TODO: change into complex_t
{
    BDDVAR n = 4;

    bool x[n];
    for(int i=0; i<(int)n; i++) x[i] = 0;
    x[0] = 1;

    MTBDD a = mtbdd_create_basis_state_double(n, x);

    VecArr_t v_arr[1 << n];

    mtbdd_to_vector_array(a, n, COLUMN_WISE_MODE, v_arr);

    print_vector_array(v_arr, n);

    test_assert(v_arr[(1 << 3)] == (double)1.0);

    return 0;
}

int
test_create_basis_state_complex()
{
    BDDVAR n = 4;

    bool x[n];
    for(int i=0; i<(int)n; i++) x[i] = 0;
    x[0] = 1;

    MTBDD a = mtbdd_create_basis_state_mpc(n, x);

    mpc_ptr v_arr[1 << n];

    mtbdd_to_vector_array_mpc(a, n, COLUMN_WISE_MODE, v_arr);

    print_vector_array_mpc(v_arr, n);

    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    test_assert(mpc_compare( (uint64_t)v_arr[(1 << 3)], (uint64_t)mpc_one));

    mpc_clear(mpc_one);

    return 0;
}

TASK_0(int, runtests)
{
    // We are not testing garbage collection
    sylvan_gc_disable();

    // Test 1
    printf("\nTesting create all zero state double.\n");
    if (test_create_all_zero_state_double()) return 1;

    printf("\nTesting create all zero state complex.\n");
    if (test_create_all_zero_state_complex()) return 1;

    // Test 2
    printf("\nTesting create basis state double.\n");
    if (test_create_basis_state_double()) return 1;

    printf("\nTesting create basis state complex.\n");
    if (test_create_basis_state_complex()) return 1;

    return 0;
}

int main()
{
    // Standard Lace initialization with 1 worker
    lace_start(1, 0);

    // Simple Sylvan initialization, also initialize BDD, MTBDD and LDD support
    sylvan_set_sizes(1LL<<20, 1LL<<20, 1LL<<16, 1LL<<16);
    sylvan_init_package(); 
    sylvan_init_mtbdd();

    printf("Mtbdd initialization complete.\n\n");

    uint32_t mpc_type = mpc_init();
    printf("Mtbdd mpc type initialization complete, mpc_type = %d.\n\n", mpc_type);

    test_assert(mpc_type == MPC_TYPE);

    int result = RUN(runtests);

    sylvan_quit();
    lace_stop();

    return result;
}









/*
//#include <stdio.h>
//#include <time.h>

//#include "qsylvan.h"

//#include <sylvan.h>
//#include <sylvan_mtbdd.h>
//#include <sylvan_edge_weights.h>

//#include "test_assert.h"

uint64_t MTBDD_WEIGHT_ZERO = 0;
uint64_t MTBDD_WEIGHT_ONE  = 1;

**
 *
 * Basis test for MTBDD method. 
 * 
/
int test_basis_state_creation_mtbdd_1()
{
    MTBDD dd0, dd1;

    bool qubits[] = {0};
    qubits[0] = 0; dd0 = mtbdd_create_basis_state(1, qubits); // |0>  function in QSimulator
    qubits[0] = 1; dd1 = mtbdd_create_basis_state(1, qubits); // |1>

    test_assert( mtbdd_countnodes(dd0) == 3  ); // Only MTBDD have more terminals, one default, one for 0 and one for 1, so 3
    test_assert( mtbdd_countnodes(dd1) == 3  );
    test_assert( mtbdd_is_unitvector(dd0, 1) );
    test_assert( mtbdd_is_unitvector(dd1, 1) );

    test_assert( mtbdd_is_ordered(dd0, 1) ); // Not needed? For sanity check!
    test_assert( mtbdd_is_ordered(dd1, 1) );

    //AMP a; // Index of edge weight, integers on the leaves (terminals)
    MTBDDMAP index; // ComplexIndex index; in case of classes // ITCN? C_index? Points only to index as integer to fraction
    qubits[0] = 0; 
    index = mtbdd_getvalue(dd0); //, qubits); 
    test_assert(index == MTBDD_WEIGHT_ONE);

    qubits[0] = 1; 
    index = mtbdd_getvalue(dd0); //, qubits); 
    test_assert(index == MTBDD_WEIGHT_ZERO);

    qubits[0] = 0; 
    index = mtbdd_getvalue(dd1); //, qubits); 
    test_assert(index == MTBDD_WEIGHT_ZERO);
    
    qubits[0] = 1; 
    index = mtbdd_getvalue(dd1); //, qubits); 
    test_assert(index == MTBDD_WEIGHT_ONE);

    printf("basis state creation mtrdd 1:     ok\n");
    return 0;
}

**
 * 
 * 
 * 
/
int test_basis_state_creation_mtbdd_2()
{
    MTBDD dd0, dd1; // Test on 3 qubits + terminal node = 4 in number

    bool qubits[] = {0, 0, 0};
    dd0 = qmdd_create_basis_state(3, qubits); // |000>

    qubits[2] = 1; 
    qubits[0] = 1; 
    dd1 = qmdd_create_basis_state(3, qubits); // |101>
    
    test_assert( mtbdd_countnodes(dd0) == 4  );
    test_assert( mtbdd_countnodes(dd1) == 4  );
    test_assert( mtbdd_is_unitvector(dd0, 3) );
    test_assert( mtbdd_is_unitvector(dd1, 3) );

    test_assert( mtbdd_is_ordered(dd0, 3) );
    test_assert( mtbdd_is_ordered(dd1, 3) );


    // for(qubits, n, dirac_str) 
    // {
    //     for(k=2; k!=0; k--)
    //     {
    //         if(dirac_str[k] == 1) 
    //             qubits[k] = 1
    //     }
    // }


    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ONE);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q2, x3); test_assert(a == AADD_ZERO);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ONE);  // |101> -> 2^2*1 + 2^1*0 + 2^0*1 = index 5, starts with 0
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);

    printf("basis state creation mtrdd 2:     ok\n");
    return 0;
}
*/