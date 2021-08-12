#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_complex_operations()
{
    complex_t ref1, ref2, ref3, ref4, val1, val2, val3, val4;
    AMP index1, index2, index3, index4;

    // comp_exact_equal / comp_approx_equal
    ref1 = comp_make((0.2+0.4), 0.0);
    ref2 = comp_make(0.6, 0.0);
    test_assert(comp_approx_equal(ref1, ref2));
    test_assert(!comp_exact_equal(ref1, ref2));

    ref1 = comp_make(1.0, 2.99999999999999855);
    ref2 = comp_make(1.0, 3.00000000000000123);
    test_assert(comp_approx_equal(ref1, ref2));
    test_assert(!comp_exact_equal(ref1, ref2));

    ref1 = comp_make((0.5+0.5), 0.0);
    ref2 = comp_make(1.0, 0.0);
    test_assert(comp_approx_equal(ref1, ref2));
    test_assert(comp_exact_equal(ref1, ref2));

    // test C_ZERO
    ref1 = comp_make(0.0, 0.0);         index1 = comp_lookup(ref1);
    ref2 = comp_make(0.0, 0.0);         index2 = comp_lookup(ref2);
    ref3 = comp_zero();                 index3 = comp_lookup(ref3);
    test_assert(index1 == index2);
    test_assert(index1 == index3);
    test_assert(index1 == C_ZERO);

    // test C_ONE
    ref1 = comp_make(1.0, 0.0);         index1 = comp_lookup(ref1);
    ref2 = comp_make(1.0, 0.0);         index2 = comp_lookup(ref2);
    test_assert(index1 == index2);
    test_assert(index1 == C_ONE);

    // Clookup, Cvalue
    ref1 = comp_make(0.5, 0.0);              index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(0.5, 0.0);              index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(1.0/flt_sqrt(2.0),0);   index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    ref4 = comp_make(1.0/flt_sqrt(2.0),0);   index4 = comp_lookup(ref4);   val4 = comp_value(index4);
    test_assert(index1 == index2);  test_assert(comp_exact_equal(val1, val2));
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    // amp_neg
    ref1 = comp_make(0.3, 4.2);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(-0.3, -4.2);       index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    index3 = amp_neg(index1);       val3 = comp_value(index3);
    index4 = amp_neg(index2);       val4 = comp_value(index4);
    test_assert(index1 == index4);  test_assert(comp_exact_equal(val1, val4));
    test_assert(index2 == index3);  test_assert(comp_exact_equal(val2, val3));

    // amp_add
    ref1 = comp_make(5.2, 1.0);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(-0.3,7.0);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(4.9, 8.0);         index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_add(index1,index2);  val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    // amp_sub
    ref1 = comp_make(1/3, 3.5);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(1/3,-1.2);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(0.0, 4.7);         index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_sub(index1,index2);  val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));
    
    // amp_mul
    ref1 = comp_make(3.0, 5.0);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(0.5, 7.0);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(-33.5, 23.5);      index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_mul(index1,index2);  val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    ref1 = comp_make(1.0/flt_sqrt(2.0),0);   index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(1.0/flt_sqrt(2.0),0);   index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(0.5, 0.0);         index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_mul(index1,index2);  val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    // amp_div
    ref1 = comp_make(1.3,-0.7);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(1.0, 0.0);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(1.3,-0.7);         index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_div(index1,index2);  val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    ref1 = comp_make(5.0, 9.0);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(-4.0,7.0);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(43./65.,-71./65.); index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    index4=amp_div(index1, index2); val4 = comp_value(index4);
    test_assert(index3 == index4);  test_assert(comp_exact_equal(val3, val4));

    if(VERBOSE) printf("complex operations:       ok\n");
    return 0;
}

int test_basis_state_creation()
{
    bool x[] = {0};

    LACE_ME;
    
    QDD q0, q1;
    x[0] = 0; q0 = qdd_create_basis_state(1, x);
    x[0] = 1; q1 = qdd_create_basis_state(1, x);
    test_assert(qdd_countnodes(q0) == 2);
    test_assert(qdd_countnodes(q1) == 2);
    test_assert(qdd_is_unitvector(q0, 1));
    test_assert(qdd_is_unitvector(q1, 1));
    test_assert(qdd_is_ordered(q0, 1));
    test_assert(qdd_is_ordered(q1, 1));

    AMP a;
    x[0] = 0; a = qdd_get_amplitude(q0, x); test_assert(a == C_ONE);
    x[0] = 1; a = qdd_get_amplitude(q0, x); test_assert(a == C_ZERO);
    x[0] = 0; a = qdd_get_amplitude(q1, x); test_assert(a == C_ZERO);
    x[0] = 1; a = qdd_get_amplitude(q1, x); test_assert(a == C_ONE);

    QDD q2, q3;
    bool x3[] = {0, 0, 0};
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q2 = qdd_create_basis_state(3, x3);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; q3 = qdd_create_basis_state(3, x3);
    test_assert(qdd_countnodes(q2) == 4);
    test_assert(qdd_countnodes(q3) == 4);
    test_assert(qdd_is_unitvector(q2, 3));
    test_assert(qdd_is_unitvector(q3, 3));
    test_assert(qdd_is_ordered(q2, 3));
    test_assert(qdd_is_ordered(q3, 3));

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ONE);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q2, x3); test_assert(a == C_ZERO);

    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == C_ZERO);

    // TODO: also test node count

    printf("basis state creation:     ok\n");
    return 0;
}

int test_vector_addition()
{ 
    QDD q0, q1, q00, q01, q10, q11, q000, q001, q010, q100;
    bool x[] = {0};
    bool x4[] = {0, 0, 0, 0};
    AMP a;
    BDDVAR nqubits;

    LACE_ME;

    // Single qubit test
    x[0] = 0; q0 = qdd_create_basis_state(1, x);
    x[0] = 1; q1 = qdd_create_basis_state(1, x);
    q00 = qdd_plus(q0, q0);
    q01 = qdd_plus(q0, q1);
    q10 = qdd_plus(q1, q0);
    q11 = qdd_plus(q1, q1);
    q000 = qdd_plus(q00, q0);
    q001 = qdd_plus(q00, q1);
    q010 = qdd_plus(q01, q0);
    q100 = qdd_plus(q10, q0);

    x[0] = 0; a = qdd_get_amplitude(q00, x); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x[0] = 1; a = qdd_get_amplitude(q00, x); test_assert(a == C_ZERO);
    x[0] = 0; a = qdd_get_amplitude(q01, x); test_assert(a == C_ONE);
    x[0] = 1; a = qdd_get_amplitude(q01, x); test_assert(a == C_ONE);
    x[0] = 0; a = qdd_get_amplitude(q10, x); test_assert(a == C_ONE);
    x[0] = 1; a = qdd_get_amplitude(q10, x); test_assert(a == C_ONE);
    x[0] = 0; a = qdd_get_amplitude(q11, x); test_assert(a == C_ZERO);
    x[0] = 1; a = qdd_get_amplitude(q11, x); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    test_assert(q01 == q10);
    test_assert(!qdd_is_unitvector(q01, 1));

    x[0] = 0; a = qdd_get_amplitude(q000, x); test_assert(a == comp_lookup(comp_make(3.0, 0.0)));
    x[0] = 1; a = qdd_get_amplitude(q000, x); test_assert(a == C_ZERO);
    x[0] = 0; a = qdd_get_amplitude(q001, x); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x[0] = 1; a = qdd_get_amplitude(q001, x); test_assert(a == C_ONE);
    x[0] = 0; a = qdd_get_amplitude(q010, x); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x[0] = 1; a = qdd_get_amplitude(q010, x); test_assert(a == C_ONE);
    x[0] = 0; a = qdd_get_amplitude(q100, x); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x[0] = 1; a = qdd_get_amplitude(q100, x); test_assert(a == C_ONE);
    test_assert(q001 == q010);
    test_assert(q001 == q100);
    test_assert(!qdd_is_unitvector(q001, 1));


    // 4 qubit test
    nqubits = 4;
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; q0 = qdd_create_basis_state(nqubits, x4);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; q1 = qdd_create_basis_state(nqubits, x4);
    q00 = qdd_plus(q0, q0);
    q01 = qdd_plus(q0, q1);
    q10 = qdd_plus(q1, q0);
    q11 = qdd_plus(q1, q1);
    q000 = qdd_plus(q00, q0);
    q001 = qdd_plus(q00, q1);
    q010 = qdd_plus(q01, q0);
    q100 = qdd_plus(q10, q0);
    test_assert(qdd_is_ordered(q000, nqubits));
    test_assert(qdd_is_ordered(q001, nqubits));
    test_assert(qdd_is_ordered(q010, nqubits));
    test_assert(qdd_is_ordered(q100, nqubits));

    // TODO: better way to assert "all others 0" --> creating function which
    // returns array x from int will help
    // q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q00, x4); test_assert(a == C_ZERO);
    test_assert(!qdd_is_unitvector(q00, 4));

    // q0 + q1 / q1 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ONE);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q01, x4); test_assert(a == C_ZERO);
    assert(q01 == q10);
    test_assert(!qdd_is_unitvector(q01, 4));

    // q1 + q1
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q11, x4); test_assert(a == C_ZERO);
    test_assert(!qdd_is_unitvector(q11, 4));

    // q0 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == comp_lookup(comp_make(3.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q000, x4); test_assert(a == C_ZERO);
    test_assert(!qdd_is_unitvector(q000, 4));

    // q0 + q0 + q1 / q0 + q1 + q0 / q1 + q0 + q0
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == comp_lookup(comp_make(2.0, 0.0)));
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 0; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ONE);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 0; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 0; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    x4[3] = 1; x4[2] = 1; x4[1] = 1; x4[0] = 1; a = qdd_get_amplitude(q001, x4); test_assert(a == C_ZERO);
    assert(q001 == q010);
    assert(q001 == q100);
    test_assert(!qdd_is_unitvector(q001, 4));

    if(VERBOSE) printf("qdd vector addition:      ok\n");
    return 0;
}

int run_qdd_tests()
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
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<11, -1, amps_backend, norm_strat);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    printf("amps backend = %d, norm strategy = %d:\n", amps_backend, norm_strat);
    int res = run_qdd_tests();

    sylvan_quit();
    lace_exit();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_stragegies; norm_strat++) {
            if (test_with(backend, norm_strat)) return 1;
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
