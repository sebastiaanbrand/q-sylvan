
// TODO: make header

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

#include "sylvan.h"
#include "test_assert.h"
#include "sylvan_int.h"

// TODO: make header

// f(x1,x2) = x1.x2 + x1'.x3 = ite(x1,x2,x3) = if x1 then x2 else x3
int
test_mtbdd_makenode_ithvar()
{
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_true ) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_false )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_true ) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_false )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_false) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_true  )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_false) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_true  )));

    return 0;
}

// Make mtbdd for f(x1,x2,x3) = x1.x2 + x3,   B->B
// Make mtbdd for f(x1,x2,x3) = x1 + x2 + x3, B->N, B->R, B->Fract, B->C
// Make mtbdd for f(x1,x2,x3) = x1.x2.x3,     B->N, B->R, B->Fract, B->C

int
test_mtbdd_makenodes_and_leafs()
{
    //
    // From: Bryant, MTBDD 1986
    //
    // f = f(x1,x2), f: Boolean -> Boolean, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                      x1
    //             x2                  x2
    //
    //          0      1            0      1
    //
    //
    //  V is a set of vertices v.
    //
    //  v is non terminal then index(v) = {1,...,n}
    //
    //  v has two children low(v) and high(v) that are vertices.
    //  if low(v) is non terminal then index(v) < index(low(v))
    //  if high(v) is non terminal then index(v) < index(high(v))
    //  v is terminal then value(v) = {0,1}
    //
    //  G has root vertex v, define recursive.
    //
    //   1/ if v is terminal, fv = 1 if value(v) = 1, fv = 0 if value(v) = 0
    //   2/ if v is non-terminal, index(v) = i, fv = xi'.f_low(v) + xi.f_high(v)
    //

    // Build up from bottom, so you can connect the returning index to the upper layer nodes.

    // Make terminals (=leafs) - layer 2 (bottom layer)
    uint32_t vartype = 1; // boolean
    uint64_t value_low = 0; 
    uint64_t value_high = 1;

    MTBDD index_leaf_low_low   = mtbdd_makeleaf(vartype, value_low);
    MTBDD index_leaf_low_high  = mtbdd_makeleaf(vartype, value_high);
    MTBDD index_leaf_high_low  = mtbdd_makeleaf(vartype, value_low);
    MTBDD index_leaf_high_high = mtbdd_makeleaf(vartype, value_high);

    // Make non terminal nodes - layer 2 - variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_low_low, index_leaf_low_high);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_high_low, index_leaf_high_high);

    // Make root node (= non terminal node) - layer 1 (top or root layer) - variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);
    printf("%ld", index_root);

    // Test leaf values
    test_assert(mtbdd_getvalue(index_leaf_low_low) == value_low);
    test_assert(mtbdd_getvalue(index_leaf_low_high) == value_high);
    test_assert(mtbdd_getvalue(index_leaf_high_low) == value_low);
    test_assert(mtbdd_getvalue(index_leaf_high_high) == value_high);

    // Test primitive functions

    // Test path length

    // Test apply method

    // Test abstract methods

    return 0;
}


int
test_mtbdd_matrix_multiplication()
{
    return 0;
}


int
test_mtbdd_apply_function()
{
    return 0;
}


// TODO: make header for test framework

TASK_0(int, runtests)
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    // Test 1
    printf("Testing mtbdd makenode and ithvar.\n");
    if (test_mtbdd_makenode_ithvar()) return 1;

    // Test 2

    return 0;
}

int main()
{
    // Standard Lace initialization with 1 worker
    lace_start(1, 0);

    // Simple Sylvan initialization, also initialize BDD, MTBDD and LDD support
    sylvan_set_sizes(1LL<<20, 1LL<<20, 1LL<<16, 1LL<<16);
    sylvan_init_package(); 
    sylvan_init_bdd();
    sylvan_init_mtbdd();   // TODO: Make argument for terminal type
    sylvan_init_ldd();

    printf("Sylvan initialization complete.\n");

    int res = RUN(runtests);

    sylvan_quit();
    lace_stop();

    return res;
}



