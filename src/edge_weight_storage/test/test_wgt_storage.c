#include <stdio.h>

#include "test_assert.h"
#include "../wgt_storage_interface.h"

bool VERBOSE = true;

int test_cmap()
{
    void *ctable = cmap_create(1<<10, 1e-14);

    uint64_t index1, index2;
    complex_t val1, val2, val3;
    int found;

    val1 = cmake(3.5, 4.7);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    for(int k=0; k<10; k++){
        found = cmap_find_or_put(ctable, &val1, &index2);
        test_assert(found == 1);
        test_assert(index2 == index2);
    }

    val1 = cmake(0.9, 2./3.);
    val2 = cmake(0.9, 2./3.);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = cmake(1.0/flt_sqrt(2.0),0); // 1/sqrt(2)
    val2 = cmake(1.0/flt_sqrt(2.0),0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *(complex_t*)cmap_get(ctable, index1);
    test_assert(flt_abs(val3.r - val1.r) < cmap_get_tolerance());
    test_assert(flt_abs(val3.i - val1.i) < cmap_get_tolerance());

    val1 = cmake(2.99999999999999855, 0.0);
    val2 = cmake(3.00000000000000123, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *(complex_t*)cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    val1 = cmake(0.0005000000000012, 0.0);
    val2 = cmake(0.0004999999999954, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);
    val3 = *(complex_t*)cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);


    // test with tolerance = 0
    cmap_free(ctable);
    ctable = cmap_create(1<<10, 0.0);

    val1 = cmake(0.9, 2./3.);
    val2 = cmake(0.9, 2./3.);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 1);
    test_assert(index1 == index2);

    val1 = cmake(2.99999999999999855, 0.0);
    val2 = cmake(3.00000000000000123, 0.0);
    found = cmap_find_or_put(ctable, &val1, &index1); test_assert(found == 0);
    found = cmap_find_or_put(ctable, &val2, &index2); test_assert(found == 0);
    test_assert(index1 != index2);
    val3 = *(complex_t*)cmap_get(ctable, index1);
    test_assert(val3.r == val1.r && val3.i == val1.i);

    cmap_free(ctable);
    if(VERBOSE) printf("cmap tests:               ok\n");
    return 0;
}


int runtests()
{
    if (test_cmap()) return 1;
    return 0;
}

int main()
{
    return runtests();
}
