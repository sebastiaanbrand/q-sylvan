#include <stdio.h>
#include <time.h>
#include <stdbool.h>

#include "util/mpreal_tree_map.h"
#include "sylvan_qdd_complex_mpreal.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_mpreal_tree_map()
{
    mpreal_tree_map_t *map = mpreal_tree_map_create(1<<10, 1e-14);

    unsigned int index1, index2;
    mpreal_complex val1, val2, val3;
    int found;

    
    val1.r = 3.5;
    val1.i = 0.0;
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = mpreal_tree_map_find_or_put(map, val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1.r = 0; val1.i = (2./3.);
    val2.r = 0; val2.i = (2./3.);
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1.r = mpfr::sqrt(0.5); val1.i = 0;
    val2.r = mpfr::sqrt(0.5); val2.i = 0;
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = mpreal_tree_map_get(map, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1.r = 0; val1.i = 2.99999999999999855;
    val2.r = 0; val2.i = 3.00000000000000123;
    test_assert(val1.i != val2.i);
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = mpreal_tree_map_get(map, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1.r = 0.00050000000000012; val1.i = 1;
    val2.r = 0.00049999999999954; val2.i = 1;
    test_assert(val1.r != val2.r);
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = mpreal_tree_map_get(map, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1.r = 14.2, val1.i = 1;
    val2.r = 14.3, val1.i = 1;
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);

    mpreal_tree_map_free(map);
    if(VERBOSE) printf("mpfr tree map tests:      ok\n");
    return 0;
}

int scope1(mpreal_tree_map_t *map)
{
    bool found;
    unsigned int index1;

    // Insert mpreal(3.5, 4.2) value into map
    mpreal_complex val1;
    val1.r = 3.5;
    val1.i = 4.2;
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpreal_tree_map_find_or_put(map, val1, &index1); test_assert(found == 1);

    return 0;
}

int scope2(mpreal_tree_map_t *map)
{
    bool found;
    unsigned int index2;

    // Can we still look up mpreal(3.5, 4.2) in the map?
    mpreal_complex val2;
    val2.r = 3.5;
    val2.i = 4.2;
    found = mpreal_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);

    return 0;
}

int test_mpfr_tree_map_scope()
{
    mpreal_tree_map_t *map = mpreal_tree_map_create(1<<10, 1e-14);
    if (scope1(map)) return 1;
    if (scope2(map)) return 1;
    mpreal_tree_map_free(map);
    if(VERBOSE) printf("mpreal scope test:        ok\n");
    return 0;
}

int test_mpreal_complex_operations()
{
    init_mpreal_amplitude_table(1<<12, 1e-14);

    // WIP: tests
    mpreal_complex ref1, ref2, ref3, ref4, val1, val2, val3, val4;
    AMP index1, index2, index3, index4;

    // comp_exact_equal / comp_approx_equal
    ref1 = mpreal_comp_make((0.2+0.4), 0.0);
    ref2 = mpreal_comp_make(0.6, 0.0);
    test_assert(mpreal_comp_approx_equal(ref1, ref2));
    test_assert(!mpreal_comp_exact_equal(ref1, ref2));

    ref1 = mpreal_comp_make(1.0, 2.99999999999999855);
    ref2 = mpreal_comp_make(1.0, 3.00000000000000123);
    test_assert(mpreal_comp_approx_equal(ref1, ref2));
    test_assert(!mpreal_comp_exact_equal(ref1, ref2));

    ref1 = mpreal_comp_make((0.5+0.5), 0.0);
    ref2 = mpreal_comp_make(1.0, 0.0);
    test_assert(mpreal_comp_approx_equal(ref1, ref2));
    test_assert(mpreal_comp_exact_equal(ref1, ref2));

    // test C_ZERO
    ref1 = mpreal_comp_make(0.0, 0.0);      index1 = mpreal_comp_lookup(ref1);
    ref2 = mpreal_comp_make(0.0, 0.0);      index2 = mpreal_comp_lookup(ref2);
    ref3 = mpreal_comp_zero();              index3 = mpreal_comp_lookup(ref3);
    test_assert(index1 == index2);
    test_assert(index1 == index3);
    test_assert(index1 == C_ZERO);

    // test C_ONE
    ref1 = mpreal_comp_make(1.0, 0.0);      index1 = mpreal_comp_lookup(ref1);
    ref2 = mpreal_comp_make(1.0, 0.0);      index2 = mpreal_comp_lookup(ref2);
    ref3 = mpreal_comp_one();               index3 = mpreal_comp_lookup(ref3);
    test_assert(index1 == index2);
    test_assert(index1 == index3);
    test_assert(index1 == C_ONE);

    // Clookup, Cvalue
    ref1 = mpreal_comp_make(0.5, 0.0);              index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(0.5, 0.0);              index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(mpreal_sqrt2(0.5),0);   index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    ref4 = mpreal_comp_make(mpreal_sqrt2(0.5),0);   index4 = mpreal_comp_lookup(ref4);   val4 = mpreal_comp_value(index4);
    test_assert(index1 == index2);  test_assert(mpreal_comp_exact_equal(val1, val2));
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    // amp_abs
    ref1 = mpreal_comp_make(-2.0, 0.0);     index1 = mpreal_comp_lookup(ref1);
    ref2 = mpreal_comp_make(2.0, 0.0);      index2 = mpreal_comp_lookup(ref2);
    index3 = mpreal_amp_abs(index1);
    test_assert(index3 == index2);
    test_assert(index3 != index1);
    ref1 = mpreal_comp_make(3.0, 4.0);      index1 = mpreal_comp_lookup(ref1);
    ref2 = mpreal_comp_make(5.0, 0.0);      index2 = mpreal_comp_lookup(ref2);
    index3 = mpreal_amp_abs(index1);
    test_assert(index3 == index2);
    test_assert(index3 != index1);

    // amp_neg
    ref1 = mpreal_comp_make(0.3, 4.2);      index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(-0.3, -4.2);    index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    index3 = mpreal_amp_neg(index1);        val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_neg(index2);        val4 = mpreal_comp_value(index4);
    test_assert(index1 == index4);  test_assert(mpreal_comp_exact_equal(val1, val4));
    test_assert(index2 == index3);  test_assert(mpreal_comp_exact_equal(val2, val3));

    // amp_add
    ref1 = mpreal_comp_make(5.2, 1.0);      index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(-0.3,7.0);      index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(4.9, 8.0);      index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_add(index1,index2); val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    // amp_sub
    ref1 = mpreal_comp_make(1/3, 3.5);      index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(1/3,-1.2);      index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(0.0, 4.7);      index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_sub(index1,index2); val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    // amp_mul
    ref1 = mpreal_comp_make(3.0, 5.0);      index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(0.5, 7.0);      index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(-33.5, 23.5);   index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_mul(index1,index2); val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    ref1 = mpreal_comp_make(1.0/mpfr::sqrt(2.0),0);     index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(1.0/mpfr::sqrt(2.0),0);     index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(0.5, 0.0);                  index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4=mpreal_amp_mul(index1,index2);               val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    // amp_div
    ref1 = mpreal_comp_make(1.3,-0.7);      index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(1.0, 0.0);      index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(1.3,-0.7);      index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_div(index1,index2); val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));

    ref1 = mpreal_comp_make(5.0, 9.0);          index1 = mpreal_comp_lookup(ref1);   val1 = mpreal_comp_value(index1);
    ref2 = mpreal_comp_make(-4.0,7.0);          index2 = mpreal_comp_lookup(ref2);   val2 = mpreal_comp_value(index2);
    ref3 = mpreal_comp_make(43./65.,-71./65.);  index3 = mpreal_comp_lookup(ref3);   val3 = mpreal_comp_value(index3);
    index4 = mpreal_amp_div(index1, index2);    val4 = mpreal_comp_value(index4);
    test_assert(index3 == index4);  test_assert(mpreal_comp_exact_equal(val3, val4));
    
    free_mpreal_amplitude_table();
    if(VERBOSE) printf("complex ops w/ mpreal:    ok\n");
    return 0;
}

int runtests()
{
    if (test_mpreal_tree_map()) return 1;
    if (test_mpfr_tree_map_scope()) return 1;
    if (test_mpreal_complex_operations()) return 1;
    return 0;
}

int main()
{
    return runtests();
}
