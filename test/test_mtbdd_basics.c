#include <stdio.h>
#include <time.h>

//#include "qsylvan.h"

#include <sylvan.h>
#include <sylvan_mtbdd.h>
#include <sylvan_edge_weights.h>

#include "test_assert.h"

uint64_t MTBDD_WEIGHT_ZERO = 0;
uint64_t MTBDD_WEIGHT_ONE  = 1;

/**
 *
 * Basis test for MTBDD method. 
 * 
*/
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

/**
 * 
 * 
 * 
*/
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

/*
    for(qubits, n, dirac_str) 
    {
        for(k=2; k!=0; k--)
        {
            if(dirac_str[k] == 1) 
                qubits[k] = 1
        }
    }
*/

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
