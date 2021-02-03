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

    // TODO: tests
    
    free_mpreal_amplitude_table();
    if(VERBOSE) printf("complex ops w/ mpreal:    WIP\n");
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
