#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_complex_operations()
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

    // test AADD_ZERO
    ref1 = cmake(0.0, 0.0);         index1 = weight_lookup(ref1);
    ref2 = cmake(0.0, 0.0);         index2 = weight_lookup(ref2);
    ref3 = czero();                 index3 = weight_lookup(ref3);
    test_assert(index1 == index2);
    test_assert(index1 == index3);
    test_assert(index1 == AADD_ZERO);

    // test AADD_ONE
    ref1 = cmake(1.0, 0.0);         index1 = weight_lookup(ref1);
    ref2 = cmake(1.0, 0.0);         index2 = weight_lookup(ref2);
    test_assert(index1 == index2);
    test_assert(index1 == AADD_ONE);

    // Clookup, Cvalue
    ref1 = cmake(0.5, 0.0);              index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(0.5, 0.0);              index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(1.0/flt_sqrt(2.0),0);   index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    ref4 = cmake(1.0/flt_sqrt(2.0),0);   index4 = weight_lookup(ref4);   weight_value(index4, &val4);
    test_assert(index1 == index2);  test_assert(weight_eq(&val1, &val2));
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_neg
    ref1 = cmake(0.3, 4.2);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(-0.3, -4.2);       index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    index3 = wgt_neg(index1);       weight_value(index3, &val3);
    index4 = wgt_neg(index2);       weight_value(index4, &val4);
    test_assert(index1 == index4);  test_assert(weight_eq(&val1, &val4));
    test_assert(index2 == index3);  test_assert(weight_eq(&val2, &val3));

    // wgt_add
    ref1 = cmake(5.2, 1.0);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(-0.3,7.0);         index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(4.9, 8.0);         index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_add(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_sub
    ref1 = cmake(1/3, 3.5);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(1/3,-1.2);         index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(0.0, 4.7);         index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_sub(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));
    
    // wgt_mul
    ref1 = cmake(3.0, 5.0);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(0.5, 7.0);         index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(-33.5, 23.5);      index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_mul(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    ref1 = cmake(1.0/flt_sqrt(2.0),0);   index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(1.0/flt_sqrt(2.0),0);   index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(0.5, 0.0);              index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_mul(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    // wgt_div
    ref1 = cmake(1.3,-0.7);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(1.0, 0.0);         index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(1.3,-0.7);         index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_div(index1,index2);  weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    ref1 = cmake(5.0, 9.0);         index1 = weight_lookup(ref1);   weight_value(index1, &val1);
    ref2 = cmake(-4.0,7.0);         index2 = weight_lookup(ref2);   weight_value(index2, &val2);
    ref3 = cmake(43./65.,-71./65.); index3 = weight_lookup(ref3);   weight_value(index3, &val3);
    index4=wgt_div(index1, index2); weight_value(index4, &val4);
    test_assert(index3 == index4);  test_assert(weight_eq(&val3, &val4));

    if(VERBOSE) printf("complex operations:       ok\n");
    return 0;
}

int test_basis_state_creation()
{
    bool x[] = {0};
    
    QMDD q0, q1;
    x[0] = 0; q0 = qmdd_create_basis_state(1, x);
    x[0] = 1; q1 = qmdd_create_basis_state(1, x);
    test_assert(aadd_countnodes(q0) == 2);
    test_assert(aadd_countnodes(q1) == 2);
    test_assert(qmdd_is_unitvector(q0, 1));
    test_assert(qmdd_is_unitvector(q1, 1));
    test_assert(aadd_is_ordered(q0, 1));
    test_assert(aadd_is_ordered(q1, 1));

    AMP a;
    x[0] = 0; a = aadd_getvalue(q0, x); test_assert(a == AADD_ONE);
    x[0] = 1; a = aadd_getvalue(q0, x); test_assert(a == AADD_ZERO);
    x[0] = 0; a = aadd_getvalue(q1, x); test_assert(a == AADD_ZERO);
    x[0] = 1; a = aadd_getvalue(q1, x); test_assert(a == AADD_ONE);

    QMDD q2, q3;
    bool x3[] = {0, 0, 0};
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q2 = qmdd_create_basis_state(3, x3);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qmdd_create_basis_state(3, x3);
    test_assert(aadd_countnodes(q2) == 4);
    test_assert(aadd_countnodes(q3) == 4);
    test_assert(qmdd_is_unitvector(q2, 3));
    test_assert(qmdd_is_unitvector(q3, 3));
    test_assert(aadd_is_ordered(q2, 3));
    test_assert(aadd_is_ordered(q3, 3));

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
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ONE);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q3, x3); test_assert(a == AADD_ZERO);

    // TODO: also test node count

    printf("basis state creation:     ok\n");
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
    q00 = aadd_plus(q0, q0);
    q01 = aadd_plus(q0, q1);
    q10 = aadd_plus(q1, q0);
    q11 = aadd_plus(q1, q1);
    q000 = aadd_plus(q00, q0);
    q001 = aadd_plus(q00, q1);
    q010 = aadd_plus(q01, q0);
    q100 = aadd_plus(q10, q0);

    x[0] = 0; a = aadd_getvalue(q00, x); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x[0] = 1; a = aadd_getvalue(q00, x); test_assert(a == AADD_ZERO);
    x[0] = 0; a = aadd_getvalue(q01, x); test_assert(a == AADD_ONE);
    x[0] = 1; a = aadd_getvalue(q01, x); test_assert(a == AADD_ONE);
    x[0] = 0; a = aadd_getvalue(q10, x); test_assert(a == AADD_ONE);
    x[0] = 1; a = aadd_getvalue(q10, x); test_assert(a == AADD_ONE);
    x[0] = 0; a = aadd_getvalue(q11, x); test_assert(a == AADD_ZERO);
    x[0] = 1; a = aadd_getvalue(q11, x); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    test_assert(q01 == q10);
    test_assert(!qmdd_is_unitvector(q01, 1));

    x[0] = 0; a = aadd_getvalue(q000, x); test_assert(a == weight_lookup(cmake(3.0, 0.0)));
    x[0] = 1; a = aadd_getvalue(q000, x); test_assert(a == AADD_ZERO);
    x[0] = 0; a = aadd_getvalue(q001, x); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x[0] = 1; a = aadd_getvalue(q001, x); test_assert(a == AADD_ONE);
    x[0] = 0; a = aadd_getvalue(q010, x); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x[0] = 1; a = aadd_getvalue(q010, x); test_assert(a == AADD_ONE);
    x[0] = 0; a = aadd_getvalue(q100, x); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x[0] = 1; a = aadd_getvalue(q100, x); test_assert(a == AADD_ONE);
    test_assert(q001 == q010);
    test_assert(q001 == q100);
    test_assert(!qmdd_is_unitvector(q001, 1));

    

    // 4 qubit test
    nqubits = 4;
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; q0 = qmdd_create_basis_state(nqubits, x4);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; q1 = qmdd_create_basis_state(nqubits, x4);
    q00 = aadd_plus(q0, q0);
    q01 = aadd_plus(q0, q1);
    q10 = aadd_plus(q1, q0);
    q11 = aadd_plus(q1, q1);
    q000 = aadd_plus(q00, q0);
    q001 = aadd_plus(q00, q1);
    q010 = aadd_plus(q01, q0);
    q100 = aadd_plus(q10, q0);
    test_assert(aadd_is_ordered(q000, nqubits));
    test_assert(aadd_is_ordered(q001, nqubits));
    test_assert(aadd_is_ordered(q010, nqubits));
    test_assert(aadd_is_ordered(q100, nqubits));

    // TODO: better way to assert "all others 0" --> creating function which
    // returns array x from int will help
    // q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q00, x4); test_assert(a == AADD_ZERO);
    test_assert(!qmdd_is_unitvector(q00, 4));

    // q0 + q1 / q1 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ONE);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q01, x4); test_assert(a == AADD_ZERO);
    test_assert(q01 == q10);
    test_assert(!qmdd_is_unitvector(q01, 4));

    // q1 + q1
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q11, x4); test_assert(a == AADD_ZERO);
    test_assert(!qmdd_is_unitvector(q11, 4));

    // q0 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == weight_lookup(cmake(3.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q000, x4); test_assert(a == AADD_ZERO);
    test_assert(!qmdd_is_unitvector(q000, 4));

    // q0 + q0 + q1 / q0 + q1 + q0 / q1 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == weight_lookup(cmake(2.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = aadd_getvalue(q001, x4); test_assert(a == AADD_ZERO);
    test_assert(q001 == q010);
    test_assert(q001 == q100);
    test_assert(!qmdd_is_unitvector(q001, 4));

    if(VERBOSE) printf("qmdd vector addition:     ok\n");
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

    return 0;
}

int test_with(int amps_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<11, -1, amps_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    printf("amps backend = %d, norm strategy = %d:\n", amps_backend, norm_strat);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
            if (test_with(backend, norm_strat)) return 1;
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
