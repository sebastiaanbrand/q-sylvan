#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include <sylvan.h>
#include <sylvan_evbdd.h>
#include <sylvan_edge_weights.h>
#include <sylvan_edge_weights_complex.h>

#include "test_assert.h"

bool VERBOSE = true;

typedef uint64_t AMP; // <- this is also defined in qsylvan.h, but that is the QSylvan layer!
typedef uint64_t QMDD;

int test_complex_operations() // No influence on DD
{
    complex_t ref1, ref2, ref3, ref4, val1, val2, val3, val4;
    AMP index1, index2, index3, index4;

    // weight_eq / weight_approx_eq
    ref1 = cmake((0.2+0.4), 0.0);
    ref2 = cmake(0.6, 0.0);
    test_assert(weight_approx_eq(&ref1, &ref2));
    test_assert(!weight_eq(&ref1, &ref2));

    ref1 = cmake(1.0, 2.99999999999999855);
    ref2 = cmake(1.0, 3.00000000000000123);
    test_assert(weight_approx_eq(&ref1, &ref2));
    test_assert(!weight_eq(&ref1, &ref2));

    ref1 = cmake((0.5+0.5), 0.0);
    ref2 = cmake(1.0, 0.0);
    test_assert(weight_approx_eq(&ref1, &ref2));
    test_assert(weight_eq(&ref1, &ref2));

    // test EVBDD_ZERO
    ref1 = cmake(0.0, 0.0);         index1 = weight_lookup(&ref1);
    ref2 = cmake(0.0, 0.0);         index2 = weight_lookup(&ref2);
    ref3 = czero();                 index3 = weight_lookup(&ref3);
    test_assert(index1 == index2);
    test_assert(index1 == index3);
    test_assert(index1 == EVBDD_ZERO);

    // test EVBDD_ONE
    ref1 = cmake(1.0, 0.0);         index1 = weight_lookup(&ref1);
    ref2 = cmake(1.0, 0.0);         index2 = weight_lookup(&ref2);
    test_assert(index1 == index2);
    test_assert(index1 == EVBDD_ONE);

    // Clookup, Cvalue
    ref1 = cmake(0.5, 0.0);              index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(0.5, 0.0);              index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(1.0/flt_sqrt(2.0),0);   index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    ref4 = cmake(1.0/flt_sqrt(2.0),0);   index4 = weight_lookup(&ref4);   weight_value(index4, &val4);
    test_assert(index1 == index2);  test_assert(weight_eq(&val1, &val2));
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_neg
    ref1 = cmake(0.3, 4.2);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(-0.3, -4.2);       index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    index3 = wgt_neg(index1);       weight_value(index3, &val3);
    index4 = wgt_neg(index2);       weight_value(index4, &val4);
    test_assert(index1 == index4);  test_assert(weight_eq(&val1, &val4));
    test_assert(index2 == index3);  test_assert(weight_eq(&val2, &val3));

    // wgt_add
    ref1 = cmake(5.2, 1.0);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(-0.3,7.0);         index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(4.9, 8.0);         index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_add(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_sub
    ref1 = cmake(1/3, 3.5);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(1/3,-1.2);         index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(0.0, 4.7);         index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_sub(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));
    
    // wgt_mul
    ref1 = cmake(3.0, 5.0);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(0.5, 7.0);         index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(-33.5, 23.5);      index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_mul(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    ref1 = cmake(1.0/flt_sqrt(2.0),0);   index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(1.0/flt_sqrt(2.0),0);   index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(0.5, 0.0);              index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_mul(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_div
    ref1 = cmake(1.3,-0.7);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(1.0, 0.0);         index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(1.3,-0.7);         index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_div(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    ref1 = cmake(5.0, 9.0);         index1 = weight_lookup(&ref1);   weight_value(index1, &val1);
    ref2 = cmake(-4.0,7.0);         index2 = weight_lookup(&ref2);   weight_value(index2, &val2);
    ref3 = cmake(43./65.,-71./65.); index3 = weight_lookup(&ref3);   weight_value(index3, &val3);
    index4=wgt_div(index1, index2); weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    if(VERBOSE) printf("complex operations:            ok\n");
    return 0;
}

int test_basis_state_creation()
{
    bool x[] = {0};
    
    QMDD q0, q1;
    x[0] = 0; q0 = qmdd_create_basis_state(1, x); // |0>
    x[0] = 1; q1 = qmdd_create_basis_state(1, x); // |1>

    test_assert(evbdd_countnodes(q0) == 2); // Only EVBDD
    test_assert(evbdd_countnodes(q1) == 2);
    test_assert(qmdd_is_unitvector(q0, 1));
    test_assert(qmdd_is_unitvector(q1, 1));
    test_assert(evbdd_is_ordered(q0, 1));
    test_assert(evbdd_is_ordered(q1, 1));
    test_assert(fabs(qmdd_get_norm(q0, 1) - 1.0) < 1e-14);
    test_assert(qmdd_get_norm(q0, 1) == 1.0);

    AMP a; // Index to edge weight
    x[0] = 0; a = evbdd_getvalue(q0, x); test_assert(a == EVBDD_ONE);
    x[0] = 1; a = evbdd_getvalue(q0, x); test_assert(a == EVBDD_ZERO);
    x[0] = 0; a = evbdd_getvalue(q1, x); test_assert(a == EVBDD_ZERO);
    x[0] = 1; a = evbdd_getvalue(q1, x); test_assert(a == EVBDD_ONE);

    QMDD q2, q3; // Test on 3 qubits + terminal node = 4 in number
    bool x3[] = {0, 0, 0};
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q2 = qmdd_create_basis_state(3, x3); // |000>
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qmdd_create_basis_state(3, x3); // |101>
    test_assert(evbdd_countnodes(q2) == 4);
    test_assert(evbdd_countnodes(q3) == 4);
    test_assert(qmdd_is_unitvector(q2, 3));
    test_assert(qmdd_is_unitvector(q3, 3));
    test_assert(evbdd_is_ordered(q2, 3));
    test_assert(evbdd_is_ordered(q3, 3));
    test_assert(fabs(qmdd_get_norm(q2, 3) - 1.0) < 1e-14);
    test_assert(qmdd_get_norm(q2, 3) == 1.0);
    test_assert(fabs(qmdd_get_norm(q3, 3) - 1.0) < 1e-14);
    test_assert(qmdd_get_norm(q3, 3) == 1.0);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ONE);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q2, x3); test_assert(a == EVBDD_ZERO);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ONE);  // |101> -> 2^2*1 + 2^1*0 + 2^0*1 = index 5, starts with 0
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = evbdd_getvalue(q3, x3); test_assert(a == EVBDD_ZERO);

    printf("basis state creation evbdd:     ok\n");
    return 0;
}

int test_vector_addition()
{ 
    QMDD q0, q1, q00, q01, q10, q11, q000, q001, q010, q100;
    bool x[] = {0};
    bool x4[] = {0, 0, 0, 0};
    AMP a;
    BDDVAR nqubits;

    // Single qubit test
    x[0] = 0; q0 = qmdd_create_basis_state(1, x);
    x[0] = 1; q1 = qmdd_create_basis_state(1, x);
    q00 = evbdd_plus(q0, q0);
    q01 = evbdd_plus(q0, q1);
    q10 = evbdd_plus(q1, q0);
    q11 = evbdd_plus(q1, q1);
    q000 = evbdd_plus(q00, q0);
    q001 = evbdd_plus(q00, q1);
    q010 = evbdd_plus(q01, q0);
    q100 = evbdd_plus(q10, q0);

    x[0] = 0; a = evbdd_getvalue(q00, x); test_assert(a == complex_lookup(2.0, 0.0));
    x[0] = 1; a = evbdd_getvalue(q00, x); test_assert(a == EVBDD_ZERO);
    x[0] = 0; a = evbdd_getvalue(q01, x); test_assert(a == EVBDD_ONE);
    x[0] = 1; a = evbdd_getvalue(q01, x); test_assert(a == EVBDD_ONE);
    x[0] = 0; a = evbdd_getvalue(q10, x); test_assert(a == EVBDD_ONE);
    x[0] = 1; a = evbdd_getvalue(q10, x); test_assert(a == EVBDD_ONE);
    x[0] = 0; a = evbdd_getvalue(q11, x); test_assert(a == EVBDD_ZERO);
    x[0] = 1; a = evbdd_getvalue(q11, x); test_assert(a == complex_lookup(2.0, 0.0));
    test_assert(q01 == q10);
    test_assert(!qmdd_is_unitvector(q01, 1));

    x[0] = 0; a = evbdd_getvalue(q000, x); test_assert(a == complex_lookup(3.0, 0.0));
    x[0] = 1; a = evbdd_getvalue(q000, x); test_assert(a == EVBDD_ZERO);
    x[0] = 0; a = evbdd_getvalue(q001, x); test_assert(a == complex_lookup(2.0, 0.0));
    x[0] = 1; a = evbdd_getvalue(q001, x); test_assert(a == EVBDD_ONE);
    x[0] = 0; a = evbdd_getvalue(q010, x); test_assert(a == complex_lookup(2.0, 0.0));
    x[0] = 1; a = evbdd_getvalue(q010, x); test_assert(a == EVBDD_ONE);
    x[0] = 0; a = evbdd_getvalue(q100, x); test_assert(a == complex_lookup(2.0, 0.0));
    x[0] = 1; a = evbdd_getvalue(q100, x); test_assert(a == EVBDD_ONE);
    test_assert(q001 == q010);
    test_assert(q001 == q100);
    test_assert(!qmdd_is_unitvector(q001, 1));



    // 4 qubit test
    nqubits = 4;
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; q0 = qmdd_create_basis_state(nqubits, x4);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; q1 = qmdd_create_basis_state(nqubits, x4);
    q00 = evbdd_plus(q0, q0);
    q01 = evbdd_plus(q0, q1);
    q10 = evbdd_plus(q1, q0);
    q11 = evbdd_plus(q1, q1);
    q000 = evbdd_plus(q00, q0);
    q001 = evbdd_plus(q00, q1);
    q010 = evbdd_plus(q01, q0);
    q100 = evbdd_plus(q10, q0);
    test_assert(evbdd_is_ordered(q000, nqubits));
    test_assert(evbdd_is_ordered(q001, nqubits));
    test_assert(evbdd_is_ordered(q010, nqubits));
    test_assert(evbdd_is_ordered(q100, nqubits));

    // TODO: better way to assert "all others 0" --> creating function which
    // returns array x from int will help
    // q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == complex_lookup(2.0, 0.0));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q00, x4); test_assert(a == EVBDD_ZERO);
    test_assert(!qmdd_is_unitvector(q00, 4));

    // q0 + q1 / q1 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ONE);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q01, x4); test_assert(a == EVBDD_ZERO);
    test_assert(q01 == q10);
    test_assert(!qmdd_is_unitvector(q01, 4));

    // q1 + q1
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == complex_lookup(2.0, 0.0));
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q11, x4); test_assert(a == EVBDD_ZERO);
    test_assert(!qmdd_is_unitvector(q11, 4));

    // q0 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == complex_lookup(3.0, 0.0));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q000, x4); test_assert(a == EVBDD_ZERO);
    test_assert(!qmdd_is_unitvector(q000, 4));

    // q0 + q0 + q1 / q0 + q1 + q0 / q1 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == complex_lookup(2.0, 0.0));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = evbdd_getvalue(q001, x4); test_assert(a == EVBDD_ZERO);
    test_assert(q001 == q010);
    test_assert(q001 == q100);
    test_assert(!qmdd_is_unitvector(q001, 4));

    if(VERBOSE) printf("evbdd vector addition:          ok\n");
    return 0;
}

int test_inner_product()
{ 
    QMDD q0, q1, q00, q01, q10, q11, q000, q001, q010, q100, q000i, q001i;
    bool x4[] = {0, 0, 0, 0};
    BDDVAR nqubits = 4;

    // vectors from test_vector_addition
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; q0 = qmdd_create_basis_state(nqubits, x4);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; q1 = qmdd_create_basis_state(nqubits, x4);
    q00 = evbdd_plus(q0, q0);    // [0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0]
    q01 = evbdd_plus(q0, q1);    // [0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0]
    q11 = evbdd_plus(q1, q1);    // [0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0]
    q000 = evbdd_plus(q00, q0);  // [0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0]
    q001 = evbdd_plus(q00, q1);  // [0 0 2 0 0 0 0 0 0 0 0 1 0 0 0 0]
    q000i = evbdd_bundle(EVBDD_TARGET(q000), wgt_mul(EVBDD_WEIGHT(q000), complex_lookup(0,1)));
                                // [0 0 3i 0 0 0 0 0 0 0 0 0 0 0 0 0]
    q001i = evbdd_bundle(EVBDD_TARGET(q001), wgt_mul(EVBDD_WEIGHT(q001), complex_lookup(0,1)));
                                // [0 0 2i 0 0 0 0 0 0 0 0 i 0 0 0 0]

    test_assert(evbdd_inner_product(q00, q00, nqubits) == complex_lookup(4.0, 0.0));
    test_assert(evbdd_inner_product(q01, q01, nqubits) == complex_lookup(2.0, 0.0));
    test_assert(evbdd_inner_product(q11, q11, nqubits) == complex_lookup(4.0, 0.0));
    test_assert(evbdd_inner_product(q000, q000, nqubits) == complex_lookup(9.0, 0.0));
    test_assert(evbdd_inner_product(q001, q001, nqubits) == complex_lookup(5.0, 0.0));
    test_assert(evbdd_inner_product(q00, q01, nqubits) == complex_lookup(2.0, 0.0));
    test_assert(evbdd_inner_product(q00, q000, nqubits) == complex_lookup(6.0, 0.0));
    test_assert(evbdd_inner_product(q01, q000, nqubits) == complex_lookup(3.0, 0.0));
    test_assert(evbdd_inner_product(q000, q001, nqubits) == complex_lookup(6.0, 0.0));
    test_assert(evbdd_inner_product(q01, q001, nqubits) == complex_lookup(3.0, 0.0));
    test_assert(evbdd_inner_product(q000i, q000i, nqubits) == complex_lookup(9.0, 0.0));
    test_assert(evbdd_inner_product(q001i, q001i, nqubits) == complex_lookup(5.0, 0.0));
    test_assert(evbdd_inner_product(q000i, q001i, nqubits) == complex_lookup(6.0, 0.0));

    if (VERBOSE) printf("evbdd inner product:            ok\n");
    return 0;
}

int run_qmdd_tests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    // basics
    if (test_complex_operations()) return 1;
    if (test_basis_state_creation()) return 1;
    if (test_vector_addition()) return 1;
    if (test_inner_product()) return 1;

    return 0;
}

int test_with(int wgt_backend, int norm_strat, int wgt_indx_bits) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<wgt_indx_bits, 1LL<<wgt_indx_bits, -1, 
                           wgt_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("wgt backend = %d, norm strat = %d, wgt indx bits = %d:\n", 
            wgt_backend, norm_strat, wgt_indx_bits);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
            if (test_with(backend, norm_strat, 11)) return 1;
            if (backend == COMP_HASHMAP) {
                // test with edge wgt index > 23 bits
                if (test_with(backend, norm_strat, 24)) return 1;
            }
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
