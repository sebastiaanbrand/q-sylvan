
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
test_mtbdd_arithmic_functions()
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
    //

    MTBDD dd1, dd2;
    MTBDD dd_plus, dd_minus, dd_times, dd_min, dd_max;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.75);
    MTBDD index_leaf_10 = mtbdd_double(0.35);
    MTBDD index_leaf_11 = mtbdd_double(0.65);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x1_low  = %ld \n", index_x1_low);
    printf("index_x1_high = %ld \n", index_x1_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    printf("index_root_node = %ld \n", index_root_node);

    dd1 = index_root_node;

    // Compute a + b
    dd_plus = mtbdd_plus(dd1, dd1);
    printf("dd_plus = %ld \n", dd_plus);
    printf("terminal 00 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_plus))));
    printf("terminal 01 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_plus))));
    printf("terminal 10 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_plus))));
    printf("terminal 11 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_plus))));

    assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_plus)))   == 0.5);
    assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_plus)))  == 1.5);
    assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_plus)))  == 0.7);
    assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_plus))) == 1.3);

    // Compute a - b
    dd_minus = mtbdd_minus(dd1, dd1);
    printf("dd_minus = %ld \n", dd_minus);
    printf("terminal 00 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_minus))));
    printf("terminal 01 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_minus))));
    printf("terminal 10 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_minus))));
    printf("terminal 11 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_minus))));

    assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_minus)))   == 0.0);
    assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_minus)))  == 0.0);
    assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_minus)))  == 0.0);
    assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_minus))) == 0.0);

    // Compute a * b
    dd_times = mtbdd_times(dd1, dd1);
    printf("dd_times = %ld \n", dd_times);
    printf("terminal 00 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_times))));
    printf("terminal 01 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_times))));
    printf("terminal 10 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_times))));
    printf("terminal 11 = %lf \n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_times))));

    assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_times)))   == 0.0625);
    assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_times)))  == 0.5625);
    //assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_times)))  == 0.1225); // Failure!?
    //assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_times))) == 0.4225); // Failure!?

    // Compute min(a, b)
    dd_min = mtbdd_min(dd1, dd1);
    assert(dd_min == dd1);         // because dd1 = dd1, no minimum of the terminals determined?

    // Define dd2 different from dd1
    index_leaf_00 = mtbdd_double(0.75);
    index_leaf_01 = mtbdd_double(0.25);
    index_leaf_10 = mtbdd_double(0.65);
    index_leaf_11 = mtbdd_double(0.35);

    // Make non-terminal nodes - middle layer, so variable x2
    index_x2 = 2;
    index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x1_low  = %ld \n", index_x1_low);
    printf("index_x1_high = %ld \n", index_x1_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);
    dd2 = index_root_node;

    // Compute min(a, b)
    dd_min = mtbdd_min(dd1, dd2);
    printf("dd_min = %ld \n", dd_min);
    printf("isleaf = %d \n", mtbddnode_isleaf(MTBDD_GETNODE(dd_min)));
    printf("type = %d \n", mtbddnode_gettype(MTBDD_GETNODE(mtbdd_getlow(mtbdd_getlow(dd_min)))));
    printf("terminal 00 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_min))));
    printf("terminal 0x = %lf \n", mtbdd_getdouble(mtbdd_getlow(dd_min)));
    printf("terminal 1x = %lf \n", mtbdd_getdouble(mtbdd_gethigh(dd_min)));
    printf("terminal  x = %lf \n", mtbdd_getdouble(dd_min));

    assert(mtbddnode_isleaf(MTBDD_GETNODE(dd_min)) == 0);
    assert(mtbdd_getdouble(mtbdd_getlow(dd_min))  == 0.25);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_min)) == 0.35);

    // Compute max(a, b)
    dd_max = mtbdd_max(dd1, dd2);
    printf("dd_max = %ld \n", dd_max);
    printf("terminal 00 = %lf \n", mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_max))));
    printf("terminal 0x = %lf \n", mtbdd_getdouble(mtbdd_getlow(dd_max)));
    printf("terminal 1x = %lf \n", mtbdd_getdouble(mtbdd_gethigh(dd_max)));
    printf("terminal  x = %lf \n", mtbdd_getdouble(dd_max));

    assert(mtbddnode_isleaf(MTBDD_GETNODE(dd_max)) == 0);
    assert(mtbdd_getdouble(mtbdd_getlow(dd_max))  == 0.75);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_max)) == 0.65);

    return 0;
}

int
test_mtbdd_abstract_plus_function_1()
{
    //
    //  Take MTBDD dd =
    //
    //                           x0
    //                      0          1
    //                           x1
    //
    //               x2                      x2
    //
    //    v1 = 0.25    v2 = 0.75   w1 = 0.35   w2 = 0.75
    //
    //  f(x1,x2) = v1 . |x1 . |x2 + v2 . |x1 . x2 + w1 . x1 . |x2 + w2 . x1 . x2
    //

    //
    //  Vector v and w addition:
    //
    //      v + w = (v1 v2)T + (w1 w2)T = (v1 + w1  v2 + w2)T
    //
    //  f(x2) = (v1 + w1).|x2.(|x1 + x1) + (v2 + w2).x2.(|x1 + x1)
    //        = (v1 + w1).|x2 + (v2 + w2).x2
    //
    //  Remove x1 by setting var_set = {1}
    //
    //  Resulted in MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //                        x2
    //
    //          v1 + w1 = 0.6    v2 + w2 = 1.5
    //

    //// Create decision diagram

    // Make f(x1,x2) as multi terminal binairy decision diagram dd
    MTBDD dd, var_set;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.75);
    MTBDD index_leaf_10 = mtbdd_double(0.35);
    MTBDD index_leaf_11 = mtbdd_double(0.75);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    dd = index_x0;

    ////  Vector v and w addition:

    // Prepare variable set to be removed from the dd with length = 1
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 1;
    //if (length_var_set > 1) var[1] = 2;
    uint32_t var_[length_var_set];

    // Test the mtbdd var_set to array and reverse function
    var_set = mtbdd_set_from_array(var, length_var_set);
    mtbdd_set_to_array(var_set, var_);

    assert(mtbdd_set_count(var_set) == length_var_set);
    assert(var[0] == var_[0]);
    if (length_var_set > 1) assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    FILE *out = fopen("..//Testing//Temporary//output_dd_plus_1.dot", "w");
    mtbdd_fprintdot(out, dd_plus);
    fclose(out);

    // Print all kinds of gets
    printf("dd_plus       = %ld\n", dd_plus);
    printf("getnumer      = %d \n", mtbdd_getnumer(dd_plus));
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_plus));
    printf("getvalue      = %ld\n", mtbdd_getvalue(dd_plus));
    printf("getlow        = %ld\n", mtbdd_getlow(dd_plus));
    printf("gethigh       = %ld\n", mtbdd_gethigh(dd_plus));
    printf("getvar        = %d \n", mtbdd_getvar(dd_plus));  // index_x2 (index_x0)

    printf("getlow(low)   00 = %lf\n", mtbdd_getdouble( mtbdd_getlow(dd_plus)));
    printf("getlow(high)  10 = %lf\n", mtbdd_getdouble( mtbdd_gethigh(dd_plus)));

    assert(mtbdd_getdouble(mtbdd_getlow(dd_plus))  == 0.6);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.5);

    return 0;

}

int
test_mtbdd_abstract_plus_function_2()
{
    //
    //  Take MTBDD dd =
    //
    //                           x0
    //                      0          1
    //                           x1
    //
    //               x2                      x2
    //
    //    v1 = 0.25    v2 = 0.75   w1 = 0.35   w2 = 0.75
    //
    //  f(x1,x2) = v1 . |x1 . |x2 + v2 . |x1 . x2 + w1 . x1 . |x2 + w2 . x1 . x2
    //

    //
    //  Vector element v and w addition:
    //
    //  f(x1) = (v1 + v2).|x1.(|x2 + x2) + (w1 + w2).x1.(|x2 + x2)
    //        = (v1 + v2).|x1 + (w1 + w2).x1
    //
    //  Remove x2 by setting var_set = {2}
    //
    //  Resulted MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //                        x1
    //
    //          v1 + v2 = 1.0    w1 + w2 = 1.1
    //

    //// Create decision diagram

    // Make f(x1,x2) as multi terminal binairy decision diagram dd
    MTBDD dd, var_set;
    //MTBDD dd_plus, dd_minus, dd_times, dd_min, dd_max;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.75);
    MTBDD index_leaf_10 = mtbdd_double(0.35);
    MTBDD index_leaf_11 = mtbdd_double(0.75);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    dd = index_x0;

    ////  Vector element v and w addition:

    // Prepare variable set to be removed from the dd with length = 1
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 2;
    //if (length_var_set > 1) var[1] = 2;
    uint32_t var_[length_var_set];
    
    // Test the mtbdd var_set to array and reverse function
    var_set = mtbdd_set_from_array(var, length_var_set);
    mtbdd_set_to_array(var_set, var_);

    assert(mtbdd_set_count(var_set) == length_var_set);
    assert(var[0] == var_[0]);
    if (length_var_set > 1) assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    FILE *out = fopen("..//Testing//Temporary//output_dd_plus_2.dot", "w");
    mtbdd_fprintdot(out, dd_plus);
    fclose(out);

    // Print all kinds of gets
    printf("dd_plus       = %ld\n", dd_plus);
    printf("getnumer      = %d \n", mtbdd_getnumer(dd_plus));
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_plus));
    printf("getvalue      = %ld\n", mtbdd_getvalue(dd_plus));
    printf("getlow        = %ld\n", mtbdd_getlow(dd_plus));
    printf("gethigh       = %ld\n", mtbdd_gethigh(dd_plus));
    printf("getvar        = %d \n", mtbdd_getvar(dd_plus));  // index_x1 (index_x0)

    printf("getlow(low)   00 = %lf\n", mtbdd_getdouble( mtbdd_getlow(dd_plus)));
    printf("getlow(high)  10 = %lf\n", mtbdd_getdouble( mtbdd_gethigh(dd_plus)));

    assert(mtbdd_getdouble(mtbdd_getlow(dd_plus))  == 1.0);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.1);

    return 0;
}

int
test_mtbdd_abstract_plus_min_max_times_function_3()
{
    //
    //  Take MTBDD dd =
    //
    //                           x0
    //                      0          1
    //                           x1
    //
    //               x2                      x2
    //
    //    v1 = 0.25    v2 = 0.75   w1 = 0.35   w2 = 0.75
    //
    //  f(x1,x2) = v1 . |x1 . |x2 + v2 . |x1 . x2 + w1 . x1 . |x2 + w2 . x1 . x2
    //

    //
    //  All vector elements v and w addition:
    //
    //  h() = f(x1) = f(x2) = v1 + v2 + w1 + w2 = v1 + w1 + v2 + w2
    //
    //  Remove x1 and x2 by setting var_set = {1,2} or {2,1}
    //
    //  Resulted MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //
    //          v1 + v2 + w1 + w2 = 2.1
    //
    //   h() = mtbdd_abstract_plus( dd, var_set = {1,2} ) 
    //

    //// Create decision diagram

    // Make f(x1,x2) as multi terminal binairy decision diagram dd
    MTBDD dd, var_set;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.75);
    MTBDD index_leaf_10 = mtbdd_double(0.35);
    MTBDD index_leaf_11 = mtbdd_double(0.75);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    dd = index_x0;

    ////  All vector elements v and w addition:

    // Prepare variable set to be removed from the dd with length = 2
    size_t length_var_set = 2;
    uint32_t var[length_var_set];
    var[0] = 1;
    if (length_var_set > 1) var[1] = 2;
    
    // Test the mtbdd var_set to array and reverse function
    uint32_t var_[length_var_set];
    var_set = mtbdd_set_from_array(var, length_var_set);
    mtbdd_set_to_array(var_set, var_);

    assert(mtbdd_set_count(var_set) == length_var_set);
    assert(var[0] == var_[0]);
    if (length_var_set > 1) assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    FILE *out = fopen("..//Testing//Temporary//output_dd_plus_3.dot", "w");
    mtbdd_fprintdot(out, dd_plus);
    fclose(out);

    // Print all kinds of gets
    printf("dd_plus       = %ld\n", dd_plus);
    printf("getnumer      = %d \n", mtbdd_getnumer(dd_plus));
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_plus));
    printf("getvalue      = %ld\n", mtbdd_getvalue(dd_plus));
    printf("getlow        = %ld\n", mtbdd_getlow(dd_plus));
    printf("gethigh       = %ld\n", mtbdd_gethigh(dd_plus));
    printf("getvar        = %d \n", mtbdd_getvar(dd_plus));  // index_x0 undefined

    assert(mtbdd_getdouble(dd_plus)  == 2.1);

    // Other operations
    MTBDD dd_min = mtbdd_abstract_min(dd, var_set);
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_min));
    assert(mtbdd_getdouble(dd_min) == 0.25);

    MTBDD dd_max = mtbdd_abstract_max(dd, var_set);
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_max));
    assert(mtbdd_getdouble(dd_max) == 0.75);

    MTBDD dd_times = mtbdd_abstract_times(dd, var_set);
    printf("getdouble     = %lf\n", mtbdd_getdouble(dd_times));
    assert(mtbdd_getdouble(dd_times) == 0.25 * 0.75 * 0.35 * 0.75 );

    return 0;
}


//
// Invent:
//
//   index_to_winning_terminal = mtbdd_path_and_index_to_minimum_leaf(dd_in, dd_out, operator) 
//
// - return index so terminal type independent
// - operator could be minimum or maximum
//
// Make wrapper to reduce operator argument in function
//
// Explanation:
//
//   mtbdd_and_abstract_plus(dd1 = A, dd2 = B, v = x2 = b)
//
// Equivalent with P(A,B), Pr(A) = sum b el B Pr(A|B=b), b = x2
//


int
test_mtbdd_and_abstract_plus_function()
{
    //
    //  Test with a matrix M[row][col] times vector v[row] multiplication
    //  
    //  M v = (m00 m01) (v0) = (m00.v0 + m01.v1) = w
    //        (m10 m11) (v1)   (m10.v0 + m11.v1)
    //
    //  Place m10 on leaf_01 and m01 on leaf_10 for M, and m01 on leaf_01 / m01 on leaf_01 for M_
    //
    //  Take MTBDD M =
    //
    //                           x0
    //                      0          1
    //                           x1
    //    col               0          1
    //                x2                      x2
    //    row    0          1          0            1
    //        
    //       m00 = 0.25  m10 = 0.35   m01 = 0.75   m11 = 0.65
    //
    //    M(x1,x2) = m00 . |x1 . |x2 + m10 . |x1 . x2 + m01 . x1 . |x2 + m11 . x1 . x2
    //
    //
    //  Take MTBDD v = 
    //
    //                           x0
    //                      0          1
    //                           x1
    //    row               0          1
    //
    //                   v0 = 3.0   v1 = 2.0
    //
    //    v(x1) = v0 . |x1 + v1 . x1
    //
    //  AND operation: 
    //
    //  w(x1,x2) = M(x1,x2) . v(x1)
    //           = m00 . v0 . |x1 . |x2 + m00 . v1 . |x1 . x1 . |x2 + m10 . v0 . |x1 . x2 + m10 . v1 . |x1 . x1 . x2 + ...
    //           = m00 . v0 . |x1 . |x2 + m10 . v0 . |x1 . x2 + m01 . v1 . x1 . |x2 + m11 . v1 . x1 . x2
    //
    //  w(x2)    = (m00 . v0 + m01 . v1) . |x2 + (m10 . v0 + m11 . v1) . x2, after elimination of x1
    //
    //  So, M v = w(x2), matrix multiplication of 2 x 2 . 2
    //
    //  Algorithm:
    //
    //    w(x2) = mtbdd_and_abstract_plus( M(x1,x2), v(x1), var_set = {1} )
    //
    //  Expected result:
    //
    //                   x0
    //                0      1
    //                   x2
    //    row         0      1
    //
    //               w0     w1
    //
    //
    //    w0 = 0.25 x 3.0 + 0.75 x 2.0
    //    w1 = 0.35 x 3.0 + 0.65 x 2.0
    //

    //// Create decision diagram for M(x1,x2)
    MTBDD M;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.35);
    MTBDD index_leaf_10 = mtbdd_double(0.75);
    MTBDD index_leaf_11 = mtbdd_double(0.65);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    M = index_x0;


    //// Create decision diagram for M_(x1,x2)
    MTBDD M_;

    // Set the terminal leafs
    index_leaf_00 = mtbdd_double(0.25);
    index_leaf_01 = mtbdd_double(0.75);
    index_leaf_10 = mtbdd_double(0.35);
    index_leaf_11 = mtbdd_double(0.65);

    printf("index_leaf_00 = %ld \n", index_leaf_00);
    printf("index_leaf_01 = %ld \n", index_leaf_01);
    printf("index_leaf_10 = %ld \n", index_leaf_10);
    printf("index_leaf_11 = %ld \n", index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    index_x2 = 2;
    index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    M_ = index_x0;


    //// Create decision diagram v(x1)
    MTBDD v;

    // Set the terminal leafs
    MTBDD index_leaf_0 = mtbdd_double(3.0);
    MTBDD index_leaf_1 = mtbdd_double(2.0);

    printf("index_leaf_0 = %ld \n", index_leaf_0);
    printf("index_leaf_1 = %ld \n", index_leaf_1);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_x0 = mtbdd_makenode(index_x1, index_leaf_0, index_leaf_1);

    printf("index_x0 = %ld \n", index_x0);

    v = index_x0;


    //// Calculate M v = w

    // Prepare variable set to be removed from the v
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 1;
    
    // Test the mtbdd var_set to array and reverse function
    MTBDD var_set = mtbdd_set_from_array(var, length_var_set);

    // Compute abstract_plus(dd, var_set)
    MTBDD w = mtbdd_and_abstract_plus(M, v, var_set);

    // Print w
    FILE *out = fopen("..//Testing//Temporary//output_and_abstract_plus_w.dot", "w");
    mtbdd_fprintdot(out, w);
    fclose(out);

    // Print all kinds of gets
    printf("w             = %ld\n", w);
    printf("getnumer      = %d \n", mtbdd_getnumer(w));
    printf("getdouble     = %lf\n", mtbdd_getdouble(w));
    printf("getvalue      = %ld\n", mtbdd_getvalue(w));
    printf("getlow        = %ld\n", mtbdd_getlow(w));
    printf("gethigh       = %ld\n", mtbdd_gethigh(w));
    printf("getvar        = %d \n", mtbdd_getvar(w));

    printf("getdouble(getlow)   0 = %lf\n", mtbdd_getdouble( mtbdd_getlow(w)  ));
    printf("getdouble(gethigh)  1 = %lf\n", mtbdd_getdouble( mtbdd_gethigh(w) ));

    double w0 = mtbdd_getdouble(mtbdd_getlow(w));
    double w1 = mtbdd_getdouble(mtbdd_gethigh(w));

    assert(mtbdd_getvar(w) == 2);
    test_assert(w0 == 0.25 * 3.0 + 0.75 * 2.0);
    test_assert(w1 == 0.35 * 3.0 + 0.65 * 2.0);


    //// Calculate M M v = M w = w_, the w has a node in x2 not x1, so use M_ because m10 <-> m01 

    // Prepare variable set to be removed from the v
    var[0] = 2;
    
    // Perpare the var_set
    var_set = mtbdd_set_from_array(var, length_var_set);
    MTBDD w_ = mtbdd_and_abstract_plus(M_, w, var_set);

    // Print w
    out = fopen("..//Testing//Temporary//output_and_abstract_plus_mw.dot", "w");
    mtbdd_fprintdot(out, w_);
    fclose(out);

    // Print all kinds of gets
    printf("w             = %ld\n", w_);
    printf("getnumer      = %d \n", mtbdd_getnumer(w_));
    printf("getdouble     = %lf\n", mtbdd_getdouble(w_));
    printf("getvalue      = %ld\n", mtbdd_getvalue(w_));
    printf("getlow        = %ld\n", mtbdd_getlow(w_));
    printf("gethigh       = %ld\n", mtbdd_gethigh(w_));
    printf("getvar        = %d \n", mtbdd_getvar(w_));

    printf("getdouble(getlow)   0 = %lf\n", mtbdd_getdouble( mtbdd_getlow(w_)  ));
    printf("getdouble(gethigh)  1 = %lf\n", mtbdd_getdouble( mtbdd_gethigh(w_) ));

    assert(mtbdd_getvar(w_) == 1);
    test_assert(mtbdd_getdouble(mtbdd_getlow(w_))  == 0.25 * w0 + 0.75 * w1);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(w_)) == 0.35 * w0 + 0.65 * w1);

    return 0;
}


int
test_mtbdd_matrix_vector_multiplication()
{
    // 
    //  M . v = w, M: 2^n x 2^n, v: 2^n x 1 (n_row x n_col)
    //

    int n = 1;

    VecArr_t v_arr[2^n];
    v_arr[0] = 1.1; v_arr[1] = -2.2;
    MTBDD v = array_vector_to_mtbdd(v_arr, n, ROW_WISE_MODE);

    // Declare and initialize two dimensional array with dynamic size
    MatArr_t **M_arr = (MatArr_t **)malloc(2^n * sizeof(MatArr_t **));
    for(int i=0; i < (2^n); i++) 
        M_arr[i] = (MatArr_t *)malloc(2^n * sizeof(MatArr_t *));

    M_arr[0][0] = -2.2; M_arr[0][1] =  1.2;
    M_arr[1][0] =  2.4; M_arr[1][1] = -1.4;
    MTBDD M = array_matrix_to_mtbdd(M_arr, n, COLUMN_WISE_MODE);

    MTBDD product = mtbdd_matvec_mult(M, v, n);

    MatArr_t w_arr[2];
    mtbdd_to_vector_array(product, n, COLUMN_WISE_MODE, w_arr);

    test_assert(w_arr[0] == M_arr[0][0] * v_arr[0] + M_arr[0][1] * v_arr[1]);
    test_assert(w_arr[1] == M_arr[1][0] * v_arr[0] + M_arr[1][1] * v_arr[1]);

    // free all mallocs
    for(int i=0; i < (2^n); i++) 
        free(M_arr[i]);

    free(M_arr);

    return 0;
}

int
test_mtbdd_matrix_matrix_multiplication()
{
    // 
    //  K . M = W, M: 2^n x 2^n, L: 2^n x 2^n, W: 2^n x 2^n
    //

    int n = 1;

    // Declare and initialize two dimensional array with dynamic size
    MatArr_t **K_arr = (MatArr_t **)malloc(2^n * sizeof(MatArr_t **));
    for(int i=0; i < (2^n); i++) 
        K_arr[i] = (MatArr_t *)malloc(2^n * sizeof(MatArr_t *));

    K_arr[0][0] =  1.1; K_arr[0][1] = 1.2; 
    K_arr[1][0] = -2.2; K_arr[1][1] = 3.3;
    MTBDD K = array_matrix_to_mtbdd(K_arr, n, ROW_WISE_MODE);

    // Declare and initialize two dimensional array with dynamic size
    MatArr_t **M_arr = (MatArr_t **)malloc(2^n * sizeof(MatArr_t **));
    for(int i=0; i < (2^n); i++) 
        M_arr[i] = (MatArr_t *)malloc(2^n * sizeof(MatArr_t *));

    M_arr[0][0] = -2.2; M_arr[0][1] =  1.2;
    M_arr[1][0] =  2.4; M_arr[1][1] = -1.4;
    MTBDD M = array_matrix_to_mtbdd(M_arr, n, COLUMN_WISE_MODE);

    MTBDD product = mtbdd_matvec_mult(K, M, n);

    // Declare and initialize two dimensional array with dynamic size
    MatArr_t **W_arr = (MatArr_t **)malloc(2^n * sizeof(MatArr_t **));
    for(int i=0; i < (2^n); i++) 
        M_arr[i] = (MatArr_t *)malloc(2^n * sizeof(MatArr_t *));

    mtbdd_to_matrix_array(product, n, COLUMN_WISE_MODE, W_arr);

    test_assert(W_arr[0][0] == M_arr[0][0] * K_arr[0][0] + M_arr[0][1] * K_arr[1][0]);
    test_assert(W_arr[0][1] == M_arr[0][0] * K_arr[0][1] + M_arr[0][1] * K_arr[1][1]);

    test_assert(W_arr[1][0] == M_arr[1][0] * K_arr[0][0] + M_arr[1][1] * K_arr[1][0]);
    test_assert(W_arr[1][1] == M_arr[1][0] * K_arr[0][1] + M_arr[1][1] * K_arr[1][1]);

    // free all mallocs
    for(int i=0; i < (2^n); i++) {
        free(K_arr[i]);
        free(M_arr[i]);
        free(W_arr[i]);
    }

    free(K_arr);
    free(M_arr);
    free(W_arr);

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
    printf("\nTesting mtbdd makenode and ithvar.\n");
    if (test_mtbdd_makenode_ithvar()) return 1;

    // Test 2
    printf("\nTesting mtbdd makeleaf, makenode, leaf type boolean integer {0,1}.\n");
    if (test_mtbdd_makenodes_and_leafs_boolean_terminals()) return 1;

    // Test 3
    printf("\nTesting mtbdd makeleaf, makenode, leaf type integer {0,1,2,3}.\n");
    if (test_mtbdd_makenodes_and_leafs_integer_terminals()) return 1;

    // Test 4
    printf("\nTesting mtbdd makeleaf, makenode, leaf type real {0.25,0.25,0.75,-0.25}.\n");
    if (test_mtbdd_makenodes_and_leafs_real_terminals()) return 1;

    // Test 5
    printf("\nTesting mtbdd arithmic functions.\n");
    if (test_mtbdd_arithmic_functions()) return 1;

    // Test 6
    printf("\nTesting mtbdd abstract arithmic functions.\n");
    if (test_mtbdd_abstract_plus_function_1()) return 1;
    if (test_mtbdd_abstract_plus_function_2()) return 1;
    if (test_mtbdd_abstract_plus_min_max_times_function_3()) return 1;

    // Test 7
    printf("\nTesting mtbdd and abstract arithmic functions.\n");
    if (test_mtbdd_and_abstract_plus_function()) return 1;

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



