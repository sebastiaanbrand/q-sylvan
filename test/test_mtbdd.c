
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
test_mtbdd_abstract_arithmic_functions()
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
    //    v1 = 0.25    v2 = 0.75   w1 = 0.35   w2 = 0.65
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
    //  Remove x2 by setting var_set = {2}
    //
    //  Resulted in MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //                        x1
    //
    //          v1 + w1 = 0.6    v2 + w2 = 1.4
    //
    //
    //  Vector element v and w addition:
    //
    //  f(x1) = (v1 + v2).|x1.(|x2 + x2) + (w1 + w2).x1.(|x2 + x2)
    //        = (v1 + v2).|x1 + (w1 + w2).x1
    //
    //  Remove x1 by setting var_set = {1}
    //
    //  Resulted MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //                        x2
    //
    //          v1 + v2 = 1.0    w1 + w2 = 1.0
    //
    //
    //  All vector elements v and w addition:
    //
    //  f() = f(x1) = f(x2) = v1 + v2 + w1 + w2 = v1 + w1 + v2 + w2
    //
    //  Remove x1 and x2 by setting var_set = {1,2} or {2,1}
    //
    //  Resulted MTBDD = mtbdd_abstract_plus(dd, var_set):
    //
    //                        x0
    //                   0          1
    //
    //          v1 + v2 + w1 + w2 = 2.0
    //

    // Make f(x1,x2) as multi terminal binairy decision diagram dd
    MTBDD dd, var_set;
    //MTBDD dd_plus, dd_minus, dd_times, dd_min, dd_max;

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
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    printf("index_x0_low  = %ld \n", index_x0_low);
    printf("index_x0_high = %ld \n", index_x0_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    printf("index_x0 = %ld \n", index_x0);

    dd = index_x0;

    // Prepare variable set to be removed from the dd with length = 1
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 1;
    if (length_var_set > 1) var[1] = 1;
    uint32_t var_[length_var_set];

    // Test the mtbdd var_set to array and reverse function
    var_set = mtbdd_set_from_array(var, length_var_set);
    mtbdd_set_to_array(var_set, var_);

    assert(mtbdd_set_count(var_set) == length_var_set);
    assert(var[0] == var_[0]);
    if (length_var_set > 1) assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);
    printf("index to result of abstract_plus = %ld \n", dd_plus);

    printf("getlow        = %ld\n", mtbdd_getlow(dd_plus));
    printf("gethigh       = %ld\n", mtbdd_gethigh(dd_plus));
    printf("getvar        = %d \n", mtbdd_getvar(dd_plus));  // index_x2 (index_x0)

/*
    printf("getlow(low)   00 = %ld\n", mtbdd_getlow(mtbdd_getlow(dd_plus)));
    printf("gethigh(low)  01 = %ld\n", mtbdd_gethigh(mtbdd_getlow(dd_plus)));
    printf("getlow(high)  10 = %ld\n", mtbdd_getlow(mtbdd_gethigh(dd_plus)));
    printf("gethigh(high) 11 = %ld\n", mtbdd_gethigh(mtbdd_gethigh(dd_plus)));
*/
    printf("getlow(low)   00 = %lf\n", mtbdd_getdouble(/*mtbdd_getlow(*/mtbdd_getlow(dd_plus)));
    //printf("gethigh(low)  01 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_plus))));
    printf("getlow(high)  10 = %lf\n", mtbdd_getdouble(/*mtbdd_getlow(*/mtbdd_gethigh(dd_plus)));
    //printf("gethigh(high) 11 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_plus))));

    //printf("getlow  = %lf\n", mtbdd_getdouble(mtbdd_getlow(dd_plus)));
    //printf("gethigh = %lf\n", mtbdd_getdouble(mtbdd_gethigh(dd_plus)));

    //assert(mtbdd_getdouble(mtbdd_getlow(dd_plus)) == 0.6);
    //assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.4);

/*
    //dd_plus = mtbdd_abstract_plus(dd1, index_x1); // Does not return in expected time, hanging ...
    //printf("index to result of abstract_plus = %ld \n", dd_plus);
    //printf("Sum of low of x1 = %lf\n", mtbdd_getdouble(mtbdd_getlow(dd_plus)));
    //printf("Sum of high of x1 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(dd_plus)));

    //assert(mtbdd_getdouble(mtbdd_getlow(dd_plus)) == 0.6);
    //assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.4);
*/
/*
    // Compute abstract_times(dd, index_of_variable)
    MTBDD dd_times = mtbdd_abstract_times(dd1, index_x2);
    printf("index to result of abstract_times = %ld \n", dd_times);
    printf("Product of low of x2 = %lf\n", mtbdd_getdouble(mtbdd_getlow(dd_times)));
    printf("Product of high of x2 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(dd_times)));

    assert(mtbdd_getdouble(mtbdd_getlow(dd_times)) == 0.0875);
    //assert(mtbdd_getdouble(mtbdd_gethigh(dd_times)) == 0.4875); // fails?

    // Compute abstract_min(dd, index_of_variable)
    MTBDD dd_min = mtbdd_abstract_min(dd1, index_x2);
    printf("index to result of abstract_min = %ld \n", dd_min);
    printf("Minimum of low of x2 = %lf\n", mtbdd_getdouble(mtbdd_getlow(dd_min)));
    printf("Minimum of high of x2 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(dd_min)));

    assert(mtbdd_getdouble(mtbdd_getlow(dd_min)) == 0.25);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_min)) == 0.65);

    // Compute abstract_max(dd, index_of_variable)
    MTBDD dd_max = mtbdd_abstract_max(dd1, index_x2);
    printf("index to result of abstract_times = %ld \n", dd_max);
    printf("Maximum of low of x2 = %lf\n", mtbdd_getdouble(mtbdd_getlow(dd_max)));
    printf("Maximum of high of x2 = %lf\n", mtbdd_getdouble(mtbdd_gethigh(dd_max)));

    assert(mtbdd_getdouble(mtbdd_getlow(dd_max)) == 0.35);
    assert(mtbdd_getdouble(mtbdd_gethigh(dd_max)) == 0.75);
*/
    return 0;

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

}


int
test_mtbdd_and_abstract_functions()
{
    //
    //  Take MTBDD mm =
    //
    //                           x0
    //                      0          1
    //                           x1
    //
    //               x2                      x2
    //
    //    m1 = 0.25    m2 = 0.75   m3 = 0.35   m4 = 0.65
    //
    //    m(x1,x2) = m1 . |x1 . |x2 + m2 . |x1 . x2 + m3 . x1 . |x2 + m4 . x1 . x2
    //
    //
    //  Take MTBDD yy = 
    //
    //                           x0
    //                      0          1
    //                           x1
    //
    //                      v1         v2
    //
    //    v(x1) = v1 . |x1 + v2 . x1
    //
    //  AND operation: 
    //
    //  h(x1,x2) = m(x1,x2) . v(x1)
    //           = m1 . v1 . |x1 . |x2 + m1 . v2 . |x1 . x1 . |x2 + m2 . v1 . |x1 . x2 + m2 . v2 . |x1 . x1 . x2 + ...
    //           = m1 . v1 . |x1 . |x2 + m2 . v1 . |x1 . x2 + m3 . v2 . x1 . |x2 + m4 . v2 . x1 . x2
    //
    //  h(x1)    = (m1 . v1 + m2 . v1) . |x1 + (m3 . v2 + m4 . v2) . x1, after elimination of x2  
    //
    //  So, M v = h(x1), matrix multiplication of 2 x 2 . 1 x 2
    //
    //  h(x1) = mtbdd_and_abstract( m(x1,x2), v(x1), var_set = {2} )
    //

 /*
    MTBDD dd1, dd2, var_set;
    //MTBDD dd_plus, dd_minus, dd_times, dd_min, dd_max;

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

    printf("index_root_node dd1 = %ld \n", index_root_node);

    dd1 = index_root_node;

    // Define dd2 different from dd1
    index_leaf_00 = mtbdd_double(0.75); // TODO change in ground terminal
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

    printf("index_root_node dd2 = %ld \n", index_root_node);

    dd2 = index_root_node;

    //var_set = index_root_node;
    //printf("%ld %ld \n", dd1, var_set); // Dummy to use dd1, var_set

    // Prepare variable set with length = 2
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 0;
    if (length_var_set > 1) var[1] = 1; 
    uint32_t var_[length_var_set];

    var_set = mtbdd_set_from_array(var, length_var_set);
    mtbdd_set_to_array(var_set, var_);

    assert(mtbdd_set_count(var_set) == length_var_set);
    assert(var[0] == var_[0]);
    if (length_var_set > 1) assert(var[1] == var_[1]);

// tests
    printf("dd1 = %ld, dd2 = %ld, var_set = %ld \n", dd1, dd2, var_set);
*/

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
    printf("\nTesting mtbdd makenode and ithvar.\n"); // TODO: print does not work?
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
    if (test_mtbdd_abstract_arithmic_functions()) return 1;

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



