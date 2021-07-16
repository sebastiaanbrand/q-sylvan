#include <stdio.h>

#include "test_assert.h"
#include "../amp_storage_interface.h"
#include "../mpfr_tree_map.h"

bool VERBOSE = true;

int test_mpfr_tree_map()
{
    mpfr_tree_map_t *map = mpfr_tree_map_create(1<<10, 1e-14);

    unsigned int index1, index2;
    mpfr_t val1, val2, sqrt_val1, sqrt_val2;
    mpfr_ptr res3;
    int found;

    mpfr_init2(val1, MPFR_PREC);
    mpfr_set_d(val1, 3.5, DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = mpfr_tree_map_find_or_put(map, val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }
    mpfr_clear(val1);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 2./3., DEFAULT_RND);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 2./3., DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    mpfr_clear(val1);
    mpfr_clear(val2);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 0.5, DEFAULT_RND);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 0.5, DEFAULT_RND);
    mpfr_init2(sqrt_val1, MPFR_PREC);
    mpfr_init2(sqrt_val2, MPFR_PREC);
    mpfr_sqrt(sqrt_val1, val1, DEFAULT_RND);
    mpfr_sqrt(sqrt_val2, val2, DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    res3 = mpfr_tree_map_get(map, index1);
    test_assert(mpfr_equal_p(res3, val1));
    mpfr_clear(val1);
    mpfr_clear(val2);
    mpfr_clear(sqrt_val1);
    mpfr_clear(sqrt_val2);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 2.999999999999999855, MPFR_PREC);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 3.000000000000000123, MPFR_PREC);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    res3 = mpfr_tree_map_get(map, index1);
    test_assert(mpfr_equal_p(res3, val1));
    mpfr_clear(val1);
    mpfr_clear(val2);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 0.0005000000000012, MPFR_PREC);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 0.0004999999999954, MPFR_PREC);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    res3 = mpfr_tree_map_get(map, index1);
    test_assert(mpfr_equal_p(res3, val1));
    mpfr_clear(val1);
    mpfr_clear(val2);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 17.000001, MPFR_PREC);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 16.999998, MPFR_PREC);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    mpfr_clear(val1);
    mpfr_clear(val2);

    // test with tolerance = 0
    mpfr_tree_map_free(map);
    map = mpfr_tree_map_create(1<<10, 0.0);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 2./3., DEFAULT_RND);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 2./3., DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    mpfr_clear(val1);
    mpfr_clear(val2);

    mpfr_init2(val1, MPFR_PREC); mpfr_set_d(val1, 0.0005000000000012, MPFR_PREC);
    mpfr_init2(val2, MPFR_PREC); mpfr_set_d(val2, 0.0004999999999954, MPFR_PREC);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    res3 = mpfr_tree_map_get(map, index1);
    test_assert(mpfr_equal_p(res3, val1));
    mpfr_clear(val1);
    mpfr_clear(val2);

    mpfr_tree_map_free(map);
    if(VERBOSE) printf("mpfr tree map tests:      ok\n");
    return 0;
}

int scope1(mpfr_tree_map_t *map)
{
    bool found;
    unsigned int index1;

    // Insert mpfr(3.5) value into map
    mpfr_t val1;
    mpfr_init2(val1, MPFR_PREC);
    mpfr_set_d(val1, 3.5, DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 0);
    found = mpfr_tree_map_find_or_put(map, val1, &index1); test_assert(found == 1);

    return 0;
}

int scope2(mpfr_tree_map_t *map)
{
    bool found;
    unsigned int index2;

    // Can we still look up mpfr(3.5) in the map?
    mpfr_t val2;
    mpfr_init2(val2, MPFR_PREC);
    mpfr_set_d(val2, 3.5, DEFAULT_RND);
    found = mpfr_tree_map_find_or_put(map, val2, &index2); test_assert(found == 1);

    return 0;
}

int test_mpfr_tree_map_scope()
{
    mpfr_tree_map_t *map = mpfr_tree_map_create(1<<10, 1e-14);
    if (scope1(map)) return 1;
    if (scope2(map)) return 1;
    mpfr_tree_map_free(map);
    if(VERBOSE) printf("mpfr scope test:          ok\n");
    return 0;
}

int runtests()
{
    if (test_mpfr_tree_map()) return 1;
    if (test_mpfr_tree_map_scope()) return 1;
    return 0;
}

int main()
{
    return runtests();
}