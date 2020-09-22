#include <stdio.h>
#include <time.h>

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_qdd_complex.h"

bool VERBOSE = true;

int test_cmap()
{
    cmap_t *ctable = cmap_create(1<<10);

    ref_t index1, index2;
    complex_t val1, val2, val3;
    int found;

    val1 = comp_make(3.5, 4.7);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = cmap_find_or_put(ctable, &val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1 = comp_make(0.9, 2./3.);
    val2 = comp_make(0.9, 2./3.);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = comp_make(1.0/sqrt(2.0),0); // 1/sqrt(2)
    val2 = comp_make(1.0/sqrt(2.0),0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *cmap_get(ctable, index1);
    test_assert(fabs(val3.r - val1.r) < TOLERANCE);
    test_assert(fabs(val3.i - val1.i) < TOLERANCE);

    val1 = comp_make(2.99999999999999855, 0.0);
    val2 = comp_make(3.00000000000000123, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1 = comp_make(0.0005000000000000012, 0.0);
    val2 = comp_make(0.0004999999999999954, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);
    
    if(VERBOSE) printf("cmap tests:               ok\n");
    cmap_free(ctable);
    return 0;
}

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
    test_assert(index1 == index2);
    test_assert(index1 == C_ZERO);

    // test C_ONE
    ref1 = comp_make(1.0, 0.0);         index1 = comp_lookup(ref1);
    ref2 = comp_make(1.0, 0.0);         index2 = comp_lookup(ref2);
    test_assert(index1 == index2);
    test_assert(index1 == C_ONE);

    // Clookup, Cvalue
    ref1 = comp_make(0.5, 0.0);         index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(0.5, 0.0);         index2 = comp_lookup(ref2);   val2 = comp_value(index2);
    ref3 = comp_make(1.0/sqrt(2.0),0);   index3 = comp_lookup(ref3);   val3 = comp_value(index3);
    ref4 = comp_make(1.0/sqrt(2.0),0);   index4 = comp_lookup(ref4);   val4 = comp_value(index4);
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

    ref1 = comp_make(1.0/sqrt(2.0),0);   index1 = comp_lookup(ref1);   val1 = comp_value(index1);
    ref2 = comp_make(1.0/sqrt(2.0),0);   index2 = comp_lookup(ref2);   val2 = comp_value(index2);
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
    
    QDD q0, q1;
    x[0] = 0; q0 = qdd_create_basis_state(1, x);
    x[0] = 1; q1 = qdd_create_basis_state(1, x);
    test_assert(qdd_countnodes(q0) == 2);
    test_assert(qdd_countnodes(q1) == 2);
    test_assert(qdd_is_unitvector(q0, 1));
    test_assert(qdd_is_unitvector(q1, 1));

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
    x4[3] = 0; x4[2] = 0; x4[1] = 1; x4[0] = 0; q0 = qdd_create_basis_state(4, x4);
    x4[3] = 1; x4[2] = 0; x4[1] = 1; x4[0] = 0; q1 = qdd_create_basis_state(4, x4);
    q00 = qdd_plus(q0, q0);
    q01 = qdd_plus(q0, q1);
    q10 = qdd_plus(q1, q0);
    q11 = qdd_plus(q1, q1);
    q000 = qdd_plus(q00, q0);
    q001 = qdd_plus(q00, q1);
    q010 = qdd_plus(q01, q0);
    q100 = qdd_plus(q10, q0);

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

int test_x_gate()
{
    QDD q0, q1, q2, q3, q4, q5;
    bool x[] = {0};
    bool x3[] = {0, 0, 0};

    LACE_ME;

    // Single qubit test
    x[0] = 0; q0 = qdd_create_basis_state(1, x);
    x[0] = 1; q1 = qdd_create_basis_state(1, x);
    x[0] = 0; q2 = qdd_create_basis_state(1, x);

    q0 = qdd_gate(q0, GATEID_X, 0); test_assert(q0 == q1);
    q0 = qdd_gate(q0, GATEID_X, 0); test_assert(q0 == q2);
    q0 = qdd_gate(q0, GATEID_X, 0); test_assert(q0 == q1);
    q0 = qdd_gate(q0, GATEID_X, 0); test_assert(q0 == q2);

    // 3 qubit test
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qdd_create_basis_state(3, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; q4 = qdd_create_basis_state(3, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; q5 = qdd_create_basis_state(3, x3);
    test_assert(qdd_countnodes(q3) == 4);
    test_assert(qdd_countnodes(q4) == 4);
    test_assert(qdd_countnodes(q5) == 4);
    
    q3 = qdd_gate(q3, GATEID_X, 1); test_assert(q3 == q4);
    q3 = qdd_gate(q3, GATEID_X, 0); test_assert(q3 == q5);
    test_assert(qdd_countnodes(q3) == 4);
    test_assert(qdd_countnodes(q4) == 4);
    test_assert(qdd_countnodes(q5) == 4);
    
    // Same 3 qubit test with sqrt(X)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; q3 = qdd_create_basis_state(3, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; q4 = qdd_create_basis_state(3, x3);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; q5 = qdd_create_basis_state(3, x3);

    q3 = qdd_gate(q3, GATEID_sqrtX, 1); test_assert(qdd_is_unitvector(q3, 3));
    q3 = qdd_gate(q3, GATEID_sqrtX, 1); test_assert(qdd_is_unitvector(q3, 3));
    test_assert(q3 == q4);

    q3 = qdd_gate(q3, GATEID_sqrtX, 0); test_assert(qdd_is_unitvector(q3, 3));
    q3 = qdd_gate(q3, GATEID_sqrtX, 0); test_assert(qdd_is_unitvector(q3, 3));
    test_assert(q3 == q5);
    
    test_assert(qdd_countnodes(q3) == 4);
    test_assert(qdd_countnodes(q4) == 4);
    test_assert(qdd_countnodes(q5) == 4);


    if(VERBOSE) printf("qdd x gates:              ok\n");
    return 0;
}

int test_h_gate()
{
    QDD q0, q1, q2, q3, q4, q5;
    bool x[] = {0};
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Single qubit test
    x[0] = 0; q0 = qdd_create_basis_state(1, x);
    x[0] = 1; q1 = qdd_create_basis_state(1, x);

    q0 = qdd_gate(q0, GATEID_H, 0);
    q1 = qdd_gate(q1, GATEID_H, 0);

    x[0] = 0; a = qdd_get_amplitude(q0, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 1; a = qdd_get_amplitude(q0, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 0; a = qdd_get_amplitude(q1, x); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x[0] = 1; a = qdd_get_amplitude(q1, x); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));


    // Two qubit test
    x2[1] = 0; x2[0] = 0; q2 = qdd_create_basis_state(2, x2); // |00>
    x2[1] = 0; x2[0] = 1; q3 = qdd_create_basis_state(2, x2); // |01>
    x2[1] = 0; x2[0] = 0; q4 = qdd_create_basis_state(2, x2); // |00>
    x2[1] = 0; x2[0] = 0; q5 = qdd_create_basis_state(2, x2); // |00>
    q2 = qdd_gate(q2, GATEID_H, 0); // q2 = |0+>
    q3 = qdd_gate(q3, GATEID_H, 0); // q3 = |0->
    q4 = qdd_gate(q4, GATEID_H, 1); // q4 = |+0>
    q5 = qdd_gate(q5, GATEID_H, 0);
    q5 = qdd_gate(q5, GATEID_H, 1); // q5 = |++>

    // q2 = |0+>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q2, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q2, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q2, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q2, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(q2) == 2);

    // q3 = |0->
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q3, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q3, x2); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q3, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q3, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(q3) == 3);

    // q4 = |+0>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q4, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q4, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q4, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q4, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(q4) == 2);

    // q5 = |++>
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q5, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(q5) == 1);

    if(VERBOSE) printf("qdd h gates:              ok\n");
    return 0;
}

int test_phase_gates()
{
    QDD q0, qZ, qS, qSS, qT, qTT, qTTTT, qTTdag, qTdagT;
    bool x2[] = {0, 0};
    AMP a;

    LACE_ME;

    // simple 2 qubit test
    x2[1] = 0; x2[0] = 0; q0 = qdd_create_basis_state(2, x2);
    q0 = qdd_gate(q0, GATEID_H, 0);
    q0 = qdd_gate(q0, GATEID_H, 1);

    qZ    = qdd_gate(q0, GATEID_Z, 0);
    qS    = qdd_gate(q0, GATEID_S, 0);
    qSS   = qdd_gate(qS, GATEID_S, 0);
    qT    = qdd_gate(q0, GATEID_T, 0);
    qTT   = qdd_gate(qT, GATEID_T, 0);
    qTTTT = qdd_gate(qTT, GATEID_T, 0);
    qTTTT = qdd_gate(qTTTT, GATEID_T, 0);
    qTTdag = qdd_gate(q0, GATEID_T, 0);
    qTTdag = qdd_gate(qTTdag, GATEID_Tdag, 0);
    qTdagT = qdd_gate(q0, GATEID_Tdag, 0);
    qTdagT = qdd_gate(qTdagT, GATEID_T, 0);

    test_assert(qZ == qSS);
    test_assert(qS == qTT);
    test_assert(qZ == qTTTT);
    test_assert(q0 == qTTdag);
    test_assert(q0 == qTdagT);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(q0) == 1);

    q0 = qdd_gate(q0, GATEID_Z, 0);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(q0) == 2);

    q0 = qdd_gate(q0, GATEID_Z, 0);
    q0 = qdd_gate(q0, GATEID_Z, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(q0) == 2);

    q0 = qdd_gate(q0, GATEID_Z, 1);
    q0 = qdd_gate(q0, GATEID_S, 0);
    q0 = qdd_gate(q0, GATEID_S, 0);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(q0) == 2);

    q0 = qdd_gate(q0, GATEID_Z, 0);
    q0 = qdd_gate(q0, GATEID_T, 1);
    q0 = qdd_gate(q0, GATEID_T, 1);
    q0 = qdd_gate(q0, GATEID_T, 1);
    q0 = qdd_gate(q0, GATEID_T, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(q0) == 2);

    q0 = qdd_gate(q0, GATEID_Z, 1);
    q0 = qdd_gate(q0, GATEID_Tdag, 1);
    q0 = qdd_gate(q0, GATEID_Tdag, 1);
    q0 = qdd_gate(q0, GATEID_Tdag, 1);
    q0 = qdd_gate(q0, GATEID_Tdag, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(q0, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(q0) == 2);


    // check R_k gates
    test_assert(gates[GATEID_Rk(0)][3] == gates[GATEID_I][3]);
    test_assert(gates[GATEID_Rk(1)][3] == gates[GATEID_Z][3]);
    test_assert(gates[GATEID_Rk(2)][3] == gates[GATEID_S][3]);
    test_assert(gates[GATEID_Rk(3)][3] == gates[GATEID_T][3]);
    test_assert(gates[GATEID_Rk_dag(0)][3] == gates[GATEID_I][3]);
    test_assert(gates[GATEID_Rk_dag(1)][3] == gates[GATEID_Z][3]);
    test_assert(gates[GATEID_Rk_dag(3)][3] == gates[GATEID_Tdag][3]);

    if(VERBOSE) printf("qdd phase gates:          ok\n");
    return 0;
}

int test_cx_gate()
{
    QDD qBell;
    bool x2[] = {0,0};
    AMP a;

    LACE_ME;

    // Test Bell state
    x2[1] = 0; x2[0] = 0; qBell = qdd_create_basis_state(2, x2);
    qBell = qdd_gate(qBell, GATEID_H, 0);
    
    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(qBell, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(qBell, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(qBell, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(qBell, x2); test_assert(a == C_ZERO);
    test_assert(qdd_countnodes(qBell) == 2);

    qBell = qdd_cgate(qBell, GATEID_X, 0, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(qBell, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(qBell, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(qBell, x2); test_assert(a == C_ZERO);
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(qBell, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    test_assert(qdd_countnodes(qBell) == 4);

    // TODO: more tests

    if(VERBOSE) printf("qdd cnot gates:           ok\n");
    return 0;
}

int test_cz_gate()
{
    QDD qGraph;
    bool x2[] = {0, 0};
    AMP a;

    LACE_ME;

    // 2 qubit graph state
    x2[1] = 0; x2[0] = 0; qGraph =qdd_create_basis_state(2, x2);
    qGraph = qdd_gate(qGraph, GATEID_H, 0);
    qGraph = qdd_gate(qGraph, GATEID_H, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    test_assert(qdd_countnodes(qGraph) == 1);

    qGraph = qdd_cgate(qGraph, GATEID_Z, 0, 1);

    x2[1] = 0; x2[0] = 0; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 0; x2[0] = 1; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 0; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(0.5, 0)));
    x2[1] = 1; x2[0] = 1; a = qdd_get_amplitude(qGraph, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
    test_assert(qdd_countnodes(qGraph) == 3);

    if(VERBOSE) printf("qdd CZ gates:             ok\n");
    return 0;
}

int test_ccz_gate()
{
    QDD q3;
    bool x3[] = {0,0,0};
    AMP a, aRef;

    LACE_ME;

    q3 = qdd_create_basis_state(3, x3);
    q3 = qdd_gate(q3, GATEID_H, 0);
    q3 = qdd_gate(q3, GATEID_H, 1);
    q3 = qdd_gate(q3, GATEID_H, 2);
    aRef = qdd_get_amplitude(q3, x3);

    x3[2]=1; x3[1]=0; x3[0]=1; q3 = qdd_all_control_phase(q3, 3, x3);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);    
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == amp_mul(aRef,comp_lookup(comp_make(-1.0,0.0))));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q3, x3); test_assert(a == aRef);

    // TODO: few more tests

    if(VERBOSE) printf("qdd all-control z gate:   ok\n");
    return 0;
}

int test_swap_circuit()
{
    QDD q;
    bool x3[] = {0,0,0};
    AMP a;

    LACE_ME;

    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_X, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    q = qdd_circuit_swap(q, 0, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    q = qdd_circuit_swap(q, 0, 1);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // TODO: more tests

    if(VERBOSE) printf("qdd swap gates:           ok\n");
    return 0;
}

int test_cswap_circuit()
{
    QDD q;
    bool x3[] = {0,0,0};
    AMP a;

    uint32_t cs[3];
    cs[0] = 0; cs[1] = QDD_INVALID_VAR; cs[2] = QDD_INVALID_VAR; // control only on q0

    LACE_ME;

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |0>, nothing should happen
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=1; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |1>, should swap q1 and q2
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); // input state: |10+> = 1/sqrt(2)(|100> + |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+>, expected output: 1/sqrt(2)(|100> + |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); 
    q = qdd_gate(q, GATEID_Z, 0); // input state: |10-> = 1/sqrt(2)(|100> - |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |->, expected output: 1/sqrt(2)(|100> - |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); 
    q = qdd_gate(q, GATEID_S, 0); // input state: |10>|+i> = 1/sqrt(2)(|100> + i|101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(0,1.0/sqrt(2.0))));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+i>, expected output: 1/sqrt(2)(|100> + i|011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(0,1.0/sqrt(2.0))));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // TODO: more tests ?

    if(VERBOSE) printf("qdd c-swap gates:         ok\n");
    return 0;
}

// measurement of q0 + sanity checks
int test_measure_random_state(QDD qdd, BDDVAR nvars)
{
    QDD qm;
    int m; double p;

    test_assert(qdd_is_unitvector(qdd, nvars));
    qm = qdd_measure_q0(qdd, nvars, &m, &p);
    test_assert(qdd_is_unitvector(qm, nvars));
    if ((fabs(p - 1.0) < TOLERANCE) || (fabs(p - 0.0) < TOLERANCE)){
        qdd = qdd_remove_global_phase(qdd); // measurement removes global phase
        test_assert(qdd_equivalent(qdd, qm, nvars, false, false));
        test_assert(qdd_equivalent(qdd, qm, nvars, true,  false));
        test_assert(qm == qdd);
    } else {
        test_assert(qdd_countnodes(qm) < qdd_countnodes(qdd));
    }

    return 0;
}

int test_measurements()
{
    QDD q, qPM;
    AMP a;
    bool x2[] = {0,0};
    bool x3[] = {0,0,0};
    int m;
    int repeat = 10;
    double prob;

    LACE_ME;
    srand(time(NULL));

    // kets labeld as |q2, q1, q0>
    for(int i=0; i < repeat; i++) {
        // |000> 
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
        x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 1);
        test_assert(prob == 0.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 1);
        qPM = qdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=m; x3[0]=0; // either |000> or |010> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 2, 3, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=m; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00->
        x3[2]=0; x3[1]=0; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        
        // |+++>, measure q0
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob); 
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=0; x3[0]=m; // either |++0> or |++1> depending on m
        q = qdd_create_basis_state(3, x3); 
        q = qdd_gate(q, GATEID_H, 1);
        q = qdd_gate(q, GATEID_H, 2);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM); 

        // |+++>, measure q1
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=m; x3[0]=0; // either |+0+> or |+1+> depending on m
        q = qdd_create_basis_state(3, x3); 
        q = qdd_gate(q, GATEID_H, 0);
        q = qdd_gate(q, GATEID_H, 2);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // [1/2, 1/2, 1/2, -1/2] = 1/2(|00> + |01> + |10> - |11>)_(q1, q0)
        x2[1]=0; x2[0]=0;
        q   = qdd_create_basis_state(2, x2);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_cgate(q,GATEID_Z, 0, 1);
        x2[1]=0; x2[0]=0; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=0; x2[0]=1; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=1; x2[0]=0; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=1; x2[0]=1; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
        qPM = qdd_measure_qubit(q, 0, 2, &m, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        if (m == 0) { // expect 1/sqrt(2)(|00> + |10>)
            x2[1]=0; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
            x2[1]=0; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=1; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
            x2[1]=1; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
        }
        if (m == 1) { // expect 1/sqrt(2)(|01> - |11>)
            x2[1]=0; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=0; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/sqrt(2.0),0)));
            x2[1]=1; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=1; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(-1.0/sqrt(2.0),0)));
        }
    }

    // Test measure all
    bool ms[3] = {0};
    int m_zer[3] = {0};
    for(int i=0; i < repeat; i++) {
        
        // |000>  
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 0);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
         x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 1);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

       // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=0; x3[0]=ms[0]; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[0] == 0) m_zer[0] += 1;

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 1);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=0; x3[1]=ms[1]; x3[0]=0; // either |000> or |010> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[1] == 0) m_zer[1] += 1;

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(fabs(prob - 0.5) < TOLERANCE);
        x3[2]=ms[2]; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[2] == 0) m_zer[2] += 1;
    }
    // test having measured |0> and |1> at least once each for |00+>, |0+0> and 
    // |+00>. (the probability that that doesn't happen is about 99.7% for 
    // repeat = 10)
    test_assert(m_zer[0] > 0);  test_assert(m_zer[0] < repeat);
    test_assert(m_zer[1] > 0);  test_assert(m_zer[1] < repeat);
    test_assert(m_zer[2] > 0);  test_assert(m_zer[2] < repeat);

    
    // TODO: more tests

    if(VERBOSE) printf("qdd measurements:         ok\n");
    return 0;
}

int test_QFT()
{
    QDD q3, q5, qref3, qref5;
    AMP a;

    // 3 qubit QFT
    bool x3[] = {0,1,1}; // little endian (q0, q1, q2)
    q3 = qdd_create_basis_state(3, x3);
    qref3 = qdd_create_basis_state(3, x3);
    q3 = qdd_circuit(q3, CIRCID_QFT, 0, 2);
    q3 = qdd_circuit(q3, CIRCID_reverse_range, 0, 2);

    // check approx equal against output from qiskit
    x3[2]=0; x3[1]=0; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(3.5355339059327384e-01,-8.6595605623549353e-17)));
    x3[2]=0; x3[1]=0; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-3.5355339059327384e-01,8.6595605623549353e-17)));
    x3[2]=0; x3[1]=1; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-1.0824450702943669e-16,-3.5355339059327384e-01)));
    x3[2]=0; x3[1]=1; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(1.0824450702943669e-16,3.5355339059327384e-01)));
    x3[2]=1; x3[1]=0; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-2.5000000000000000e-01,2.5000000000000017e-01)));
    x3[2]=1; x3[1]=0; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(2.5000000000000000e-01,-2.5000000000000017e-01)));
    x3[2]=1; x3[1]=1; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(2.5000000000000017e-01,2.5000000000000000e-01)));
    x3[2]=1; x3[1]=1; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-2.5000000000000017e-01,-2.5000000000000000e-01)));
    test_assert(qdd_is_unitvector(q3, 3));

    // inverse QFT
    q3 = qdd_circuit(q3, CIRCID_reverse_range, 0, 2);
    q3 = qdd_circuit(q3, CIRCID_QFT_inv, 0, 2);
    test_assert(qdd_equivalent(q3, qref3, 3, false, false));
    test_assert(qdd_equivalent(q3, qref3, 3, true, false));
    test_assert(q3 == qref3);


    // 5 qubit QFT
    bool x5[] = {0,1,1,0,1};
    q5 = qdd_create_basis_state(5, x5);
    qref5 = qdd_create_basis_state(5, x5);
    q5 = qdd_circuit(q5, CIRCID_QFT, 0, 4);
    q5 = qdd_circuit(q5, CIRCID_reverse_range, 0, 4);

    // check approx equal against output from qiskit
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7677669529663692e-01,-6.4946704217662027e-17)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7677669529663692e-01,6.4946704217662027e-17)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(7.5771154920605696e-17,1.7677669529663692e-01)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-7.5771154920605696e-17,-1.7677669529663692e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.2500000000000011e-01,-1.2499999999999999e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.2500000000000011e-01,1.2499999999999999e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.2499999999999999e-01,-1.2500000000000011e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.2499999999999999e-01,1.2500000000000011e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(6.7649512518274585e-02,-1.6332037060954713e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-6.7649512518274585e-02,1.6332037060954713e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.6332037060954713e-01,6.7649512518274571e-02)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.6332037060954713e-01,-6.7649512518274571e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.6332037060954710e-01,6.7649512518274696e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.6332037060954710e-01,-6.7649512518274696e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-6.7649512518274710e-02,-1.6332037060954710e-01)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(6.7649512518274710e-02,1.6332037060954710e-01)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.4698445030241986e-01,9.8211869798387877e-02)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.4698445030241986e-01,-9.8211869798387877e-02)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-9.8211869798387877e-02,-1.4698445030241986e-01)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(9.8211869798387877e-02,1.4698445030241986e-01)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7337998066526852e-01,3.4487422410367806e-02)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7337998066526852e-01,-3.4487422410367806e-02)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-3.4487422410367799e-02,1.7337998066526852e-01)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(3.4487422410367799e-02,-1.7337998066526852e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(3.4487422410367972e-02,1.7337998066526850e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-3.4487422410367972e-02,-1.7337998066526850e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7337998066526850e-01,3.4487422410367986e-02)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7337998066526850e-01,-3.4487422410367986e-02)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(9.8211869798387752e-02,-1.4698445030241994e-01)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-9.8211869798387752e-02,1.4698445030241994e-01)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.4698445030241994e-01,9.8211869798387752e-02)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.4698445030241994e-01,-9.8211869798387752e-02)));
    test_assert(qdd_is_unitvector(q5, 5));

    // inverse QFT
    q5 = qdd_circuit(q5, CIRCID_reverse_range, 0, 4);
    q5 = qdd_circuit(q5, CIRCID_QFT_inv, 0, 4);
    test_assert(qdd_equivalent(q5, qref5, 5, false, false));
    test_assert(qdd_equivalent(q5, qref5, 5, true, false));
    test_assert(q5 == qref5);
    
    
    if(VERBOSE) printf("qdd QFT:                  ok\n");
    return 0;
}


int test_shor()
{   
    QDD q, qref;
    bool x3[] = {0,0,0};
    bool as[] = {1,0,1};
    AMP a;

    LACE_ME;

    // <Test qdd_phi_add>
    // Test inversion
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0);
    q = qdd_gate(q, GATEID_H, 1);
    q = qdd_gate(q, GATEID_H, 2);
    qref = q;
    q = qdd_phi_add(q, 0, 2, as);
    q = qdd_phi_add_inv(q, 0, 2, as);
    test_assert(qdd_equivalent(q, qref, 3, false, false));
    test_assert(qdd_equivalent(q, qref, 3, true, false));
    test_assert(q == qref);

    // Test addition in Fourier space (be mindful about endianness here!)
    // 2 + 1 (no carry)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 0; as[2] = 0;          // a = 100 = 1 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |3> = |011> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 2 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |4> = |100> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 2 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 0;          // x = 010 = 2 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 3 + 2 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 0; as[1] = 1; as[2] = 0;          // a = 010 = 2 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |5> = |101> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // 3 + 3 (carry, should go the the left)
    x3[0] = 0; x3[1] = 1; x3[2] = 1;          // x = 011 = 3 (MSB first)
    as[0] = 1; as[1] = 1; as[2] = 0;          // a = 110 = 3 (LSB first)
    q = qdd_create_basis_state(3, x3);        // create state |x>
    q = qdd_circuit(q, CIRCID_QFT, 0, 2);     // QFT|x> = |phi(x)>
    q = qdd_phi_add(q, 0, 2, as);             // addition in Fourier space gives |phi(x+a)>
    q = qdd_circuit(q, CIRCID_QFT_inv, 0, 2); // expected out = |6> = |110> (MSB first)
    x3[0]=0; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=0; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=0; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[0]=1; x3[1]=1; x3[2]=0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[0]=1; x3[1]=1; x3[2]=1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    // </Test qdd_phi_add>

    if(VERBOSE) printf("qdd fourier addition:     ok\n");


    // test Shor
    srand(time(NULL));
    uint64_t N, factor, counter, nqubits;

    // 15 = 3 x 5 (11 qubits)
    N = 15;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = run_shor(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qdd %ld-qubit Shor:        ok (found factor %ld of %ld with %ld tries)\n", nqubits, factor, N, counter);

    // 35 = 5 x 7 (15 qubits)
    N = 35;
    counter = 0;
    factor = 0;
    nqubits = ceil(log2(N))*2 + 3;
    while (!factor) {
        factor = run_shor(N, 0, false);
        counter++;
    }
    if(VERBOSE) printf("qdd %ld-qubit Shor:        ok (found factor %ld of %ld with %ld tries)\n", nqubits, factor, N, counter);


    return 0;
}


int test_5qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    int n_qubits = 5;
    bool x5[] = {0,0,0,0,0};
    FILE *fp;
    
    LACE_ME;

    // 5 qubit state
    qref = qdd_create_basis_state(n_qubits, x5);
    q    = qdd_create_basis_state(n_qubits, x5);

    // 32 gates 
    q = qdd_cgate(q, GATEID_Z, 1, 2);       q = qdd_gate(q, GATEID_X, 2);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 3, 4);       q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 3);           q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_cgate(q, GATEID_X, 1, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_cgate(q, GATEID_X, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 0, 1);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 2, 3);       q = qdd_gate(q, GATEID_Z, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 3, 4);       q = qdd_cgate(q, GATEID_Z, 0, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_Z, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 3);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_gate(q, GATEID_H, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 3);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 1, 2);       q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, n_qubits)) return 1;

    // inverse
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_H, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_X, 0, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_gate(q, GATEID_Z, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 2);       q = qdd_cgate(q, GATEID_X, 3, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 0);           q = qdd_cgate(q, GATEID_Z, 2, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 1);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 4);       q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 2);       q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_gate(q, GATEID_H, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_cgate(q, GATEID_Z, 3, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qdd_is_unitvector(q, 5));

    test_assert(qdd_equivalent(q, qref, n_qubits, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, n_qubits, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    fp = fopen("5_qubit_res.dot", "w");
    qdd_fprintdot(fp, q, true);
    fclose(fp);

    if(VERBOSE) printf("qdd 5 qubit circuit:      ok (%ld nodes)\n", node_count);
    return 0;
}

int test_10qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    bool x10[] = {0,0,0,0,0,0,0,0,0,0};
    
    LACE_ME;

    // 10 qubit state
    qref = qdd_create_basis_state(10, x10);
    q    = qdd_create_basis_state(10, x10);

    // 30 random* Clifford gates            *chosen by a fair dice roll
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_X, 6, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 3, 5);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 3, 8);       q = qdd_cgate(q, GATEID_Z, 3, 6);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 0, 7);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 1, 9);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 2);       q = qdd_gate(q, GATEID_X, 2);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_X, 0, 4);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 8);       q = qdd_cgate(q, GATEID_X, 6, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 9);       q = qdd_gate(q, GATEID_X, 9);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 4, 9);       q = qdd_cgate(q, GATEID_Z, 2, 7);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 7, 8);       q = qdd_gate(q, GATEID_X, 7);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 2);           q = qdd_gate(q, GATEID_Z, 7);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 10));
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, 10)) return 1;

    // inverse
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_X, 6);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 7);           q = qdd_gate(q, GATEID_Z, 2);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 7);           q = qdd_cgate(q, GATEID_Z, 7, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 2, 7);       q = qdd_cgate(q, GATEID_X, 4, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_cgate(q, GATEID_X, 0, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 6, 9);       q = qdd_cgate(q, GATEID_X, 0, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 4);       q = qdd_cgate(q, GATEID_X, 5, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_X, 0, 2);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 1, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_gate(q, GATEID_Z, 3);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 3, 6);       q = qdd_cgate(q, GATEID_X, 3, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_H, 1);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 6, 9);       q = qdd_gate(q, GATEID_X, 6);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 10));

    test_assert(qdd_equivalent(q, qref, 10, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, 10, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qdd 10 qubit circuit:     ok (%ld nodes)\n", node_count);
    return 0;
}

int test_20qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    bool x20[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    LACE_ME;

    // 20 qubit state
    qref = qdd_create_basis_state(20, x20);
    q    = qdd_create_basis_state(20, x20);

    // 100 gates
    q = qdd_cgate(q, GATEID_Z, 4, 18);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_X, 1, 12);      q = qdd_gate(q, GATEID_Z, 4);
    q = qdd_cgate(q, GATEID_Z, 9, 10);      q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 1, 16);      q = qdd_cgate(q, GATEID_Z, 13, 16);
    q = qdd_cgate(q, GATEID_Z, 7, 11);      q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_cgate(q, GATEID_Z, 1, 4);       q = qdd_cgate(q, GATEID_X, 6, 16);
    q = qdd_cgate(q, GATEID_X, 3, 18);      q = qdd_cgate(q, GATEID_X, 2, 15);      q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_gate(q, GATEID_Z, 6);
    q = qdd_cgate(q, GATEID_X, 3, 6);       q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_gate(q, GATEID_Z, 18);
    q = qdd_cgate(q, GATEID_Z, 14, 15);     q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_gate(q, GATEID_H, 8);           q = qdd_gate(q, GATEID_X, 9);
    q = qdd_gate(q, GATEID_X, 8);           q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_X, 17);          q = qdd_gate(q, GATEID_Z, 11);
    q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 4);           q = qdd_gate(q, GATEID_X, 18);
    q = qdd_cgate(q, GATEID_X, 4, 10);      q = qdd_gate(q, GATEID_X, 15);          q = qdd_cgate(q, GATEID_Z, 16, 18);     q = qdd_cgate(q, GATEID_Z, 0, 15);
    q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 16);          q = qdd_cgate(q, GATEID_X, 7, 18);
    q = qdd_gate(q, GATEID_X, 16);          q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_X, 9, 10);      q = qdd_gate(q, GATEID_X, 6);
    q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 11);          q = qdd_cgate(q, GATEID_Z, 4, 5);       q = qdd_gate(q, GATEID_X, 1);
    q = qdd_cgate(q, GATEID_Z, 2, 19);      q = qdd_cgate(q, GATEID_X, 8, 9);       q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_Z, 11, 16);
    q = qdd_cgate(q, GATEID_X, 13, 19);     q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_gate(q, GATEID_X, 6);           q = qdd_gate(q, GATEID_X, 15);
    q = qdd_gate(q, GATEID_Z, 0);           q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_gate(q, GATEID_Z, 8);
    q = qdd_cgate(q, GATEID_X, 12, 14);     q = qdd_cgate(q, GATEID_Z, 2, 18);      q = qdd_cgate(q, GATEID_X, 12, 15);     q = qdd_gate(q, GATEID_X, 9);
    q = qdd_gate(q, GATEID_Z, 12);          q = qdd_gate(q, GATEID_X, 3);           q = qdd_gate(q, GATEID_X, 0);           q = qdd_cgate(q, GATEID_X, 1, 4);
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_gate(q, GATEID_X, 19);          q = qdd_gate(q, GATEID_X, 5);           q = qdd_cgate(q, GATEID_Z, 2, 16);
    q = qdd_gate(q, GATEID_X, 4);           q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_gate(q, GATEID_Z, 12);
    q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_gate(q, GATEID_Z, 13);          q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_Z, 10);
    q = qdd_gate(q, GATEID_X, 4);           q = qdd_gate(q, GATEID_Z, 16);          q = qdd_cgate(q, GATEID_Z, 4, 17);      q = qdd_gate(q, GATEID_Z, 7);
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_Z, 6, 7);       q = qdd_cgate(q, GATEID_X, 12, 19);     q = qdd_gate(q, GATEID_Z, 15);
    q = qdd_cgate(q, GATEID_X, 5, 11);      q = qdd_cgate(q, GATEID_X, 9, 17);      q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 11, 18);
    q = qdd_cgate(q, GATEID_Z, 5, 15);      q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_cgate(q, GATEID_X, 1, 6);       q = qdd_cgate(q, GATEID_X, 8, 16);
    q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_cgate(q, GATEID_Z, 3, 18);      q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_Z, 14, 18);
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, 20)) return 1;

    // inverse
    q = qdd_cgate(q, GATEID_Z, 14, 18);     q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_Z, 3, 18);      q = qdd_cgate(q, GATEID_X, 5, 19);
    q = qdd_cgate(q, GATEID_X, 8, 16);      q = qdd_cgate(q, GATEID_X, 1, 6);       q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_cgate(q, GATEID_Z, 5, 15);
    q = qdd_cgate(q, GATEID_X, 11, 18);     q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 9, 17);      q = qdd_cgate(q, GATEID_X, 5, 11);
    q = qdd_gate(q, GATEID_Z, 15);          q = qdd_cgate(q, GATEID_X, 12, 19);     q = qdd_cgate(q, GATEID_Z, 6, 7);       q = qdd_gate(q, GATEID_H, 4);
    q = qdd_gate(q, GATEID_Z, 7);           q = qdd_cgate(q, GATEID_Z, 4, 17);      q = qdd_gate(q, GATEID_Z, 16);          q = qdd_gate(q, GATEID_X, 4);
    q = qdd_gate(q, GATEID_Z, 10);          q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_Z, 13);          q = qdd_cgate(q, GATEID_X, 9, 11);
    q = qdd_gate(q, GATEID_Z, 12);          q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_gate(q, GATEID_X, 4);
    q = qdd_cgate(q, GATEID_Z, 2, 16);      q = qdd_gate(q, GATEID_X, 5);           q = qdd_gate(q, GATEID_X, 19);          q = qdd_gate(q, GATEID_H, 1);
    q = qdd_cgate(q, GATEID_X, 1, 4);       q = qdd_gate(q, GATEID_X, 0);           q = qdd_gate(q, GATEID_X, 3);           q = qdd_gate(q, GATEID_Z, 12);
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_cgate(q, GATEID_X, 12, 15);     q = qdd_cgate(q, GATEID_Z, 2, 18);      q = qdd_cgate(q, GATEID_X, 12, 14);
    q = qdd_gate(q, GATEID_Z, 8);           q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_gate(q, GATEID_Z, 0);
    q = qdd_gate(q, GATEID_X, 15);          q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_cgate(q, GATEID_X, 13, 19);
    q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_X, 8, 9);       q = qdd_cgate(q, GATEID_Z, 2, 19);
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_cgate(q, GATEID_Z, 4, 5);       q = qdd_gate(q, GATEID_Z, 11);          q = qdd_gate(q, GATEID_X, 18);
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_X, 9, 10);      q = qdd_gate(q, GATEID_X, 2);           q = qdd_gate(q, GATEID_X, 16);
    q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_Z, 16);          q = qdd_gate(q, GATEID_X, 18);          q = qdd_cgate(q, GATEID_X, 7, 10);
    q = qdd_cgate(q, GATEID_Z, 0, 15);      q = qdd_cgate(q, GATEID_Z, 16, 18);     q = qdd_gate(q, GATEID_X, 15);          q = qdd_cgate(q, GATEID_X, 4, 10);
    q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 4);           q = qdd_gate(q, GATEID_X, 18);          q = qdd_cgate(q, GATEID_X, 12, 16);
    q = qdd_gate(q, GATEID_Z, 11);          q = qdd_gate(q, GATEID_X, 17);          q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_X, 8);
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_gate(q, GATEID_H, 8);           q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_Z, 14, 15);
    q = qdd_gate(q, GATEID_Z, 18);          q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_X, 3, 6);
    q = qdd_gate(q, GATEID_Z, 6);           q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_cgate(q, GATEID_X, 2, 15);      q = qdd_cgate(q, GATEID_X, 3, 18);
    q = qdd_cgate(q, GATEID_X, 6, 16);      q = qdd_cgate(q, GATEID_Z, 1, 4);       q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_cgate(q, GATEID_Z, 7, 11);
    q = qdd_cgate(q, GATEID_Z, 13, 16);     q = qdd_cgate(q, GATEID_Z, 1, 16);      q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 9, 10);
    q = qdd_gate(q, GATEID_Z, 4);           q = qdd_cgate(q, GATEID_X, 1, 12);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_Z, 4, 18);

    test_assert(qdd_equivalent(q, qref, 20, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, 20, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qdd 20 qubit circuit:     ok (%ld nodes)\n", node_count);
    return 0;
}



int runtests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    if (test_cmap()) return 1;
    if (test_complex_operations()) return 1;
    if (test_basis_state_creation()) return 1;
    if (test_vector_addition()) return 1;
    if (test_x_gate()) return 1;
    if (test_h_gate()) return 1;
    if (test_phase_gates()) return 1;
    if (test_cx_gate()) return 1;
    if (test_cz_gate()) return 1;
    if (test_ccz_gate()) return 1;
    if (test_swap_circuit()) return 1;
    if (test_cswap_circuit()) return 1;
    if (test_measurements()) return 1;
    if (test_5qubit_circuit()) return 1;
    if (test_10qubit_circuit()) return 1;
    if (test_20qubit_circuit()) return 1;
    if (test_QFT()) return 1;
    //if (test_grover()) return 1;
    if (test_shor()) return 1;

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
    sylvan_init_qdd(1LL<<19);

    int res = runtests();

    sylvan_quit();
    lace_exit();

    return res;
}
