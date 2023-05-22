
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
test_mtbdd_makenodes_and_leafs_boolean_terminals()
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
    //  v is terminal then value(v) = {0,1} in case of a boolean terminal.
    //
    //  G has root vertex v, define recursive.
    //
    //   1/ if v is terminal, fv = 1 if value(v) = 1, fv = 0 if value(v) = 0
    //   2/ if v is non-terminal, index(v) = i, fv = xi'.f_low(v) + xi.f_high(v)
    //
    // Built test-MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                     x1
    //
    //                  0      1
    //

    // Make terminals (=leafs) - layer 3 (bottom layer)
    uint32_t terminal_type = 0;  // terminal has integer type

    // Set the terminal leafs
    uint64_t value_low_00  = 0;
    uint64_t value_high_01 = 1;
    uint64_t value_low_10  = 0;
    uint64_t value_high_11 = 1;

    MTBDD index_leaf_00 = mtbdd_makeleaf(terminal_type, value_low_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(terminal_type, value_high_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(terminal_type, value_low_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(terminal_type, value_high_11);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Identical terminals must have the same index
    test_assert(index_leaf_00 == index_leaf_10);
    test_assert(index_leaf_01 == index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x1_low  = %ld \n", index_x1_low);
    printf("index_x1_high = %ld \n", index_x1_high);

    // The indices of x1 should be identical
    test_assert(index_x1_low == index_x1_high);

    // Make non-terminal nodes - root layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    printf("index_root_node = %ld \n", index_root_node);

    // The index of root should be the indices of x1
    test_assert(index_root_node == index_x1_low);
    test_assert(index_root_node == index_x1_high);

    return 0;
}

int
test_mtbdd_makenodes_and_leafs_integer_terminals()
{
    //
    // From: Bryant, MTBDD 1986
    //
    // f = f(x1,x2), f: Boolean -> Boolean, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                      x1
    //             x2                  x2
    //
    //          0      1            3      4
    //
    //
    // Built test-MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                      x1
    //             x2                 x2
    //
    //          0      1           3      4
    //

    // Make terminals (=leafs)
    uint32_t terminal_type = 0;  // terminal has integer type

    // Set the terminal leafs
    uint64_t value_low_00  = 0;
    uint64_t value_high_01 = 1;
    uint64_t value_low_10  = 2;
    uint64_t value_high_11 = 3;

    MTBDD index_leaf_00 = mtbdd_makeleaf(terminal_type, value_low_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(terminal_type, value_high_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(terminal_type, value_low_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(terminal_type, value_high_11);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Different terminals should have different indices
    test_assert(index_leaf_00 != index_leaf_10); 
    test_assert(index_leaf_01 != index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x1_low  = %ld \n", index_x1_low);
    printf("index_x1_high = %ld \n", index_x1_high);

    // The indices of x1 should be different
    test_assert(index_x1_low != index_x1_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    printf("index_root_node = %ld \n", index_root_node);

    // The index of root should be the indices of x1
    test_assert(index_root_node != index_x1_low);
    test_assert(index_root_node != index_x1_high);

    // Check the leaf values
    test_assert(mtbdd_getvalue(index_leaf_00) == value_low_00);
    test_assert(mtbdd_getvalue(index_leaf_01) == value_high_01);
    test_assert(mtbdd_getvalue(index_leaf_10) == value_low_10);
    test_assert(mtbdd_getvalue(index_leaf_11) == value_high_11);

    // Check of node type being leaf of terminals
    test_assert(mtbdd_isleaf(index_leaf_00) == (int)1); // (int)1 == true, TODO: make boolean!
    test_assert(mtbdd_isleaf(index_leaf_01) == (int)1);
    test_assert(mtbdd_isleaf(index_leaf_10) == (int)1);
    test_assert(mtbdd_isleaf(index_leaf_11) == (int)1);

    // Check of node type being non-terminal
    test_assert(mtbdd_isleaf(index_x1) == (int)0);
    test_assert(mtbdd_isleaf(index_x2) == (int)0);
    test_assert(mtbdd_isleaf(index_root_node) == (int)0);

    return 0;
}

int
test_mtbdd_makenodes_and_leafs_real_terminals()
{
    //
    // From: Bryant, MTBDD 1986
    //
    // f = f(x1,x2), f: Boolean -> Boolean, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                      x1
    //             x2                  x2
    //
    //        0.25      0.25     0.75     -0.25
    //
    //
    // Built test-MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                      x1
    //                              x2
    //
    //            0.25          0.75      -0.25
    //

    // Make terminals (=leafs)
    //uint32_t terminal_type = 2;  // terminal has real type

    // Set the terminal leafs
    double value_low_00  =  0.25;
    double value_high_01 =  0.25;
    double value_low_10  =  0.75;
    double value_high_11 = -0.25;

    MTBDD index_leaf_00 = mtbdd_double(/*terminal_type,*/ value_low_00);
    MTBDD index_leaf_01 = mtbdd_double(/*terminal_type,*/ value_high_01);
    MTBDD index_leaf_10 = mtbdd_double(/*terminal_type,*/ value_low_10);
    MTBDD index_leaf_11 = mtbdd_double(/*terminal_type,*/ value_high_11);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Different terminals should have different indices
    test_assert(index_leaf_00 == index_leaf_01); 
    test_assert(index_leaf_01 != index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x1_low  = %ld \n", index_x1_low);
    printf("index_x1_high = %ld \n", index_x1_high);

    // The indices of x1 should be different
    test_assert(index_x1_low != index_x1_high);
    test_assert(index_x1_low == index_leaf_00);
    test_assert(index_x1_low == index_leaf_01);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    printf("index_root_node = %ld \n", index_root_node);

    // The index of root should be the indices of x1
    test_assert(index_root_node != index_x1_low);
    test_assert(index_root_node != index_x1_high);

    // Check the leaf values
    test_assert(mtbdd_getdouble(index_leaf_00) == value_low_00);
    test_assert(mtbdd_getdouble(index_leaf_01) == value_high_01);
    test_assert(mtbdd_getdouble(index_leaf_10) == value_low_10);
    test_assert(mtbdd_getdouble(index_leaf_11) == value_high_11);

    // Check of node type being leaf of terminals
    test_assert(mtbdd_isleaf(index_leaf_00) == (int)1); // (int)1 == true, TODO: make boolean!
    test_assert(mtbdd_isleaf(index_leaf_01) == (int)1);
    test_assert(mtbdd_isleaf(index_leaf_10) == (int)1);
    test_assert(mtbdd_isleaf(index_leaf_11) == (int)1);

    printf("index_x1_low    = %d \n", mtbdd_isleaf(index_x1_low));
    printf("index_x1_high   = %d \n", mtbdd_isleaf(index_x1_high));
    printf("index_root_node = %d \n", mtbdd_isleaf(index_root_node));

    // Check of node type being non-terminal
    test_assert(mtbdd_isleaf(index_x1_low) == (int)1); // This should be identical with index_leaf_0x
    test_assert(mtbdd_isleaf(index_x1_high) == (int)0);
    test_assert(mtbdd_isleaf(index_root_node) == (int)0);

    return 0;
}

// <- above all okay
/*
    // - test getlow/high of a node

    printf("test: getlow/high of a node \n");
    printf("getlow(index_root_node = %ld) = %ld == index_x1_low = %ld \n", index_root_node, mtbdd_getlow(index_root_node), index_x1_low);
    
    test_assert(mtbdd_getlow(index_root_node) == index_x1_low);
    test_assert(mtbdd_gethigh(index_root_node) == index_x1_high);

// <- above fails ...

    // - test get type/value of terminals
    test_assert(mtbdd_gettype(index_leaf_00) == terminal_type);  // leaf, should be integer 
    test_assert(mtbdd_getvalue(index_leaf_01) == 1);             // leaf, should be 1

// <- above works

    // - test get type/value of non terminals

    printf("gettype(index_x1_low)  = %d \n", mtbdd_gettype(index_x1_low));
    printf("getvalue(index_x1_low) = %ld \n", mtbdd_getvalue(index_x1_low));

    test_assert(mtbdd_gettype(index_x1_low) == 1);  // not a leaf, 
    test_assert(mtbdd_getvalue(index_x1_low) == 0); // not a leaf, TODO: what to return if no value?

    // - test get var of non terminals
    test_assert(mtbdd_getvar(index_x1_low) == index_x2);
    test_assert(mtbdd_getvar(index_x1_high) == index_x2);
    test_assert(mtbdd_getvar(index_root_node) == index_x1);

// <- above fails ...


    // Test path length //

    // Test apply method

    // Test abstract methods

    return 0;
}
*/


int
test_mtbdd_apply_function()
{
    //
    // MTBDD dd1 =
    //
    //              x1
    //
    //       x2             x2
    //
    //  0.25    0.75    0.35    0.65
    //
    //
    // MTBDD dd2 = 
    //
    //              x1
    //
    //        x2            x2
    //
    //  0.75    0.25    0.65    0.35
    //
    //

    MTBDD dd1, dd2;

    

    mtbdd_apply(dd1, dd2, plus);

    return 0;
}

int
test_mtbdd_matrix_multiplication()
{
    // 
    //  K . L = M
    //
    //  K = (1.0  2.0)   L = (1.0  0.5)   M = (1.0 x 1.0 + 2.0 x 0.5  1.0 x 1.0 + 2.0 x 1.0)
    //      (2.0  1.0)       (0.5  1.0)       (2.0 x 1.0 + 1.0 x 1.0  2.0 x 0.5 + 1.0 x 1.0)
    //
    //  M = (2.0  3.0)
    //      (2.0  2.0)
    //

    return 0;
}

int
test_mtbdd_matrix_kronecker_multiplication()
{
    //
    //  K (x) L = M
    //
    //  K = (1.0  2.0)   L = (1.0  0.5)   M = (1.0 x L  2.0 x L)
    //      (2.0  1.0)       (0.5  1.0)       (2.0 x L  1.0 x L)
    //
    //  M = (1.0 0.5 2.0 1.0)
    //      (0.5 1.0 1.0 2.0)
    //      (2.0 1.0 1.0 0.5)
    //      (1.0 2.0 0.5 1.0)
    //

    return 0;
}


// TODO: make header for test framework

TASK_0(int, runtests)
{
    // We are not testing garbage collection
    sylvan_gc_disable();

    // Test 1
    printf("Testing mtbdd makenode and ithvar.\n"); // TODO: print does not work?
    if (test_mtbdd_makenode_ithvar()) return 1;

    // Test 2
    printf("Testing mtbdd makeleaf, makenode, leaf type boolean integer {0,1}.\n");
    if (test_mtbdd_makenodes_and_leafs_boolean_terminals()) return 1;

    // Test 3
    printf("Testing mtbdd makeleaf, makenode, leaf type integer {0,1,2,3}.\n");
    if (test_mtbdd_makenodes_and_leafs_integer_terminals()) return 1;

    // Test 4
    printf("Testing mtbdd makeleaf, makenode, leaf type real {0.25,0.25,0.75,-0.25}.\n");
    if (test_mtbdd_makenodes_and_leafs_real_terminals()) return 1;

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



