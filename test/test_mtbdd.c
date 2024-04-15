
/*
 * Copyright 2023 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

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
#include "sylvan_mpc.h"

// Testing mtbdd makenode and ithvar.

/**
 * f(x1,x2) = x1.x2 + x1'.x3 = ite(x1,x2,x3) = if x1 then x2 else x3
 */ 
int
test_mtbdd_makenode_ithvar()
{
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_true ) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_false )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_true ) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_false )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_false) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_true  )));
    test_assert(mtbdd_makenode(mtbdd_ithvar(1), sylvan_false, sylvan_false) == sylvan_not(mtbdd_makenode(mtbdd_ithvar(1), sylvan_true,  sylvan_true  )));

    return 0;
}

// Testing mtbdd makeleaf, makenode, leaf type boolean, integer, double, complex.

int
test_mtbdd_makenodes_and_leafs_boolean()
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
    //  G has root vertex v, defined recursive.
    //
    //   1/ if v is terminal, fv = 1 if value(v) = 1, fv = 0 if value(v) = 0
    //   2/ if v is non-terminal, index(v) = i, fv = xi'.f_low(v) + xi.f_high(v)
    //
    // Built MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                     x1
    //
    //                  0      1
    //

    // Make terminals (=leafs) - layer 3 (bottom layer)
    uint32_t terminal_type = 0;  // terminal has boolean type

    // Set the terminal leafs
    uint64_t value_low_00  = 0;
    uint64_t value_high_01 = 1;
    uint64_t value_low_10  = 0;
    uint64_t value_high_11 = 1;

    MTBDD index_leaf_00 = mtbdd_makeleaf(terminal_type, value_low_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(terminal_type, value_high_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(terminal_type, value_low_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(terminal_type, value_high_11);

    if(false) {
        printf("index_leaf_00 = %"PRIu64" \n", index_leaf_00);
        printf("index_leaf_01 = %"PRIu64" \n", index_leaf_01);
        printf("index_leaf_10 = %"PRIu64" \n", index_leaf_10);
        printf("index_leaf_11 = %"PRIu64" \n", index_leaf_11);
    }

    // Identical terminals must have the same index
    test_assert(index_leaf_00 == index_leaf_10);
    test_assert(index_leaf_01 == index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    if(false) {
        printf("index_x1_low  = %"PRIu64" \n", index_x1_low);
        printf("index_x1_high = %"PRIu64" \n", index_x1_high);
    }

    // The indices of x1 should be identical
    test_assert(index_x1_low == index_x1_high);

    // Make non-terminal nodes - root layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    if(false) {
        printf("index_root_node = %"PRIu64" \n", index_root_node);
    }

    // The index of root should be the indices of x1
    test_assert(index_root_node == index_x1_low);
    test_assert(index_root_node == index_x1_high);

    return 0;
}

int
test_mtbdd_makenodes_and_leafs_integer()
{
    //
    // f = f(x1,x2), f: Boolean -> Integer, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                      x1
    //             x2                  x2
    //
    //          0      1            3      4
    //
    //
    // Built MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
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

    if(false) {
        printf("index_leaf_00 = %"PRIu64" \n", index_leaf_00);
        printf("index_leaf_01 = %"PRIu64" \n", index_leaf_01);
        printf("index_leaf_10 = %"PRIu64" \n", index_leaf_10);
        printf("index_leaf_11 = %"PRIu64" \n", index_leaf_11);
    }

    // Different terminals should have different indices
    test_assert(index_leaf_00 != index_leaf_10); 
    test_assert(index_leaf_01 != index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    if(false) {
        printf("index_x1_low  = %"PRIu64" \n", index_x1_low);
        printf("index_x1_high = %"PRIu64" \n", index_x1_high);
    }

    // The indices of x1 should be different
    test_assert(index_x1_low != index_x1_high);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    if(false) {
        printf("index_root_node = %"PRIu64" \n", index_root_node);
    }

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
test_mtbdd_makenodes_and_leafs_double()
{
    //
    // f = f(x1,x2), f: Boolean -> Double, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                      x1
    //             x2                  x2
    //
    //        0.25      0.25     0.75     -0.25
    //
    //
    // Built MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                      x1
    //                              x2
    //
    //            0.25          0.75      -0.25
    //

    // Set the terminal leafs
    double value_low_00  =  0.25;
    double value_high_01 =  0.25;
    double value_low_10  =  0.75;
    double value_high_11 = -0.25;

    MTBDD index_leaf_00 = mtbdd_double(value_low_00);
    MTBDD index_leaf_01 = mtbdd_double(value_high_01);
    MTBDD index_leaf_10 = mtbdd_double(value_low_10);
    MTBDD index_leaf_11 = mtbdd_double(value_high_11);

    // Different terminals should have different indices
    test_assert(index_leaf_00 == index_leaf_01); 
    test_assert(index_leaf_01 != index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // The indices of x1 should be different
    test_assert(index_x1_low != index_x1_high);
    test_assert(index_x1_low == index_leaf_00);
    test_assert(index_x1_low == index_leaf_01);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

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

    // Check of node type being non-terminal
    test_assert(mtbdd_isleaf(index_x1_low) == (int)1); // This should be identical with index_leaf_0x
    test_assert(mtbdd_isleaf(index_x1_high) == (int)0);
    test_assert(mtbdd_isleaf(index_root_node) == (int)0);

    return 0;
}

int
test_mtbdd_makenodes_and_leafs_complex()
{
    //
    // f = f(x1,x2), f: Boolean -> Complex, f = x1 <op> f(x2) + !x1 <op> f(x2)
    //
    //                                    x1
    //                      x2                         x2
    //
    //        0.25 + 0.25i      0.25 + 0.25i     0.75     -0.25
    //
    //
    // Built MTBDD up from bottom, so you can connect the returning index to the upper layer nodes.
    //
    // The above diagram will be reduced while building up, and should result in:
    //
    //                      x1
    //                                      x2
    //
    //           0.25 + 0.25i         0.75      -0.25
    //

    // Set the terminal leafs
    mpc_t value_00;

    //mpc_assign(value_00, 0.25, 0.25);

    mpc_init2(value_00, MPC_PRECISION);
    mpc_set_d_d(value_00, 0.25, 0.25, MPC_ROUNDING);

    mpc_t value_01;
    mpc_init2(value_01, MPC_PRECISION);
    mpc_set_d_d(value_01, 0.25, 0.25, MPC_ROUNDING);

    mpc_t value_10;
    mpc_init2(value_10, MPC_PRECISION);
    mpc_set_d_d(value_10, 0.75, 0.0, MPC_ROUNDING);

    mpc_t value_11;
    mpc_init2(value_11, MPC_PRECISION);
    mpc_set_d_d(value_11, -0.25, 0.0, MPC_ROUNDING);

    if(false) {
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, value_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, value_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, value_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, value_11, MPC_ROUNDING);
        putchar('\n');
    }

    MTBDD index_leaf_00 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_00);
    MTBDD index_leaf_00_cp = mtbdd_makeleaf(MPC_TYPE, (size_t)value_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_11);

    // Different terminals should have different indices
    test_assert(index_leaf_00 == index_leaf_00_cp);
    test_assert(index_leaf_00 == index_leaf_01); 
    test_assert(index_leaf_01 != index_leaf_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // The indices of x1 should be different
    test_assert(index_x1_low != index_x1_high);
    test_assert(index_x1_low == index_leaf_00);
    test_assert(index_x1_low == index_leaf_01);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    // The index of root should be the indices of x1
    test_assert(index_root_node != index_x1_low);
    test_assert(index_root_node != index_x1_high);

    // Check the leaf values
    test_assert( mpc_compare( mtbdd_getvalue(index_leaf_00), (uint64_t)value_00) );
    test_assert( mpc_compare( mtbdd_getvalue(index_leaf_01), (uint64_t)value_01) );
    test_assert( mpc_compare( mtbdd_getvalue(index_leaf_10), (uint64_t)value_10) );
    test_assert( mpc_compare( mtbdd_getvalue(index_leaf_11), (uint64_t)value_11) );

    // Check of node type being leaf of terminals
    test_assert(mtbdd_isleaf(index_leaf_00) == (int)true);
    test_assert(mtbdd_isleaf(index_leaf_01) == (int)true);
    test_assert(mtbdd_isleaf(index_leaf_10) == (int)true);
    test_assert(mtbdd_isleaf(index_leaf_11) == (int)true);

    // Check of node type being non-terminal
    test_assert(mtbdd_isleaf(index_x1_low) == (int)1);       // This should be identical with index_leaf_0x
    test_assert(mtbdd_isleaf(index_x1_high) == (int)0);
    test_assert(mtbdd_isleaf(index_root_node) == (int)0);

    // Clear the allocated memory
    mpc_clear(value_00);
    mpc_clear(value_01);
    mpc_clear(value_10);
    mpc_clear(value_11);

    return 0;
}

// Testing mtbdd arithmic functions double / complex

int
test_mtbdd_arithmic_functions_double()
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
    // Test plus, min, times, min, max on dd1 dd1 and dd1 dd2
    //

    MTBDD dd1, dd2;
    MTBDD dd_plus, dd_minus, dd_times, dd_min, dd_max;

    // Set the terminal leafs
    double value_00 = 0.25;
    double value_01 = 0.75;
    double value_10 = 0.35;
    double value_11 = 0.65;

    MTBDD index_leaf_00 = mtbdd_double(value_00);
    MTBDD index_leaf_01 = mtbdd_double(value_01);
    MTBDD index_leaf_10 = mtbdd_double(value_10);
    MTBDD index_leaf_11 = mtbdd_double(value_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    dd1 = index_root_node;

    // Compute a + b
    dd_plus = mtbdd_plus(dd1, dd1);
    
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_plus)))   == value_00+value_00);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_plus)))  == value_01+value_01);
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_plus)))  == value_10+value_10);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_plus))) == value_11+value_11);

    // Compute a - b
    dd_minus = mtbdd_minus(dd1, dd1);
    
    test_assert(mtbdd_getdouble(mtbdd_getlow( mtbdd_getlow(dd_minus) )) == value_00-value_00);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_minus) )) == value_01-value_01);
    test_assert(mtbdd_getdouble(mtbdd_getlow( mtbdd_gethigh(dd_minus))) == value_10-value_10);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_minus))) == value_11-value_11);

    // Compute a * b
    dd_times = mtbdd_times(dd1, dd1);
    
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(dd_times)))   == value_00*value_00);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_getlow(dd_times)))  == value_01*value_01);
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(dd_times)))  == value_10*value_10);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(mtbdd_gethigh(dd_times))) == value_11*value_11);

    // Compute min(a, b)
    dd_min = mtbdd_min(dd1, dd1);
    test_assert(dd_min == dd1);

    // Define dd2 different from dd1
    double value_00_ = 0.75;
    double value_01_ = 0.25;
    double value_10_ = 0.65;
    double value_11_ = 0.35;

    index_leaf_00 = mtbdd_double(value_00_);
    index_leaf_01 = mtbdd_double(value_01_);
    index_leaf_10 = mtbdd_double(value_10_);
    index_leaf_11 = mtbdd_double(value_11_);

    // Make non-terminal nodes - middle layer, so variable x2
    index_x2 = 2;
    index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);
    dd2 = index_root_node;

    // Compute min(a, b)
    dd_min = mtbdd_min(dd1, dd2);

    test_assert(mtbddnode_isleaf(MTBDD_GETNODE(dd_min)) == 0);
    test_assert(mtbdd_getdouble(mtbdd_getlow(dd_min))  == 0.25);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(dd_min)) == 0.35);

    // Compute max(a, b)
    dd_max = mtbdd_max(dd1, dd2);
    
    test_assert(mtbddnode_isleaf(MTBDD_GETNODE(dd_max)) == 0);
    test_assert(mtbdd_getdouble(mtbdd_getlow(dd_max))  == 0.75);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(dd_max)) == 0.65);

    return 0;
}

int
test_mtbdd_arithmic_plus_sub_times_functions_complex()
{
    //
    // MTBDD dd1 =
    //
    //                         x1
    //
    //            x2                       x2
    //
    //  0.25+1.25i  0.75+1.75i    0.35+1.35i  0.65+1.65i
    //


    MTBDD dd1;
    MTBDD dd_plus, dd_times, dd_minus;

    // Set the terminal leafs
    mpc_t value_00;
    mpc_assign(value_00, 0.25, 1.25);
    mpc_t value_01;
    mpc_assign(value_01, 0.75, 1.75);
    mpc_t value_10;
    mpc_assign(value_10, 0.35, 1.35);
    mpc_t value_11;
    mpc_assign(value_11, 0.65, 1.65);

    MTBDD index_leaf_00 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    dd1 = index_root_node;

    // Compute a + b
    //dd_plus = mpc_plus(dd1, dd1);
    dd_plus = mtbdd_plus(dd1, dd1);

    mpc_t add_00;
    mpc_init2(add_00, MPC_PRECISION);
    mpc_add(add_00, value_00, value_00, MPC_ROUNDING);
    mpc_t add_01;
    mpc_init2(add_01, MPC_PRECISION);
    mpc_add(add_01, value_01, value_01, MPC_ROUNDING);
    mpc_t add_10;
    mpc_init2(add_10, MPC_PRECISION);
    mpc_add(add_10, value_10, value_10, MPC_ROUNDING);
    mpc_t add_11;
    mpc_init2(add_11, MPC_PRECISION);
    mpc_add(add_11, value_11, value_11, MPC_ROUNDING);

    if(false) {
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)add_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)add_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)add_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)add_11, MPC_ROUNDING);
        putchar('\n');
    }

    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow(mtbdd_getlow(  dd_plus))), (uint64_t)add_00));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_getlow( dd_plus))), (uint64_t)add_01));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow(mtbdd_gethigh( dd_plus))), (uint64_t)add_10));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_gethigh(dd_plus))), (uint64_t)add_11));

    mpc_clear(add_00);
    mpc_clear(add_01);
    mpc_clear(add_10);
    mpc_clear(add_11);

    // Compute a . b
    //dd_times = mpc_times(dd1, dd1);
    dd_times = mtbdd_times(dd1, dd1);

    mpc_t times_00;
    mpc_init2(times_00, MPC_PRECISION);
    mpc_mul(times_00, value_00, value_00, MPC_ROUNDING);
    mpc_t times_01;
    mpc_init2(times_01, MPC_PRECISION);
    mpc_mul(times_01, value_01, value_01, MPC_ROUNDING);
    mpc_t times_10;
    mpc_init2(times_10, MPC_PRECISION);
    mpc_mul(times_10, value_10, value_10, MPC_ROUNDING);
    mpc_t times_11;
    mpc_init2(times_11, MPC_PRECISION);
    mpc_mul(times_11, value_11, value_11, MPC_ROUNDING);

    if(false) {
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)times_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)times_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)times_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)times_11, MPC_ROUNDING);
        putchar('\n');
    }

    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow(mtbdd_getlow(  dd_times))), (uint64_t)times_00));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_getlow( dd_times))), (uint64_t)times_01));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow(mtbdd_gethigh( dd_times))), (uint64_t)times_10));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_gethigh(dd_times))), (uint64_t)times_11));

    mpc_clear(times_00);
    mpc_clear(times_01);
    mpc_clear(times_10);
    mpc_clear(times_11);

    // Compute a - b
    dd_minus = mpc_minus(dd1, dd1);

    mpc_t sub_00;
    mpc_init2(sub_00, MPC_PRECISION);
    mpc_sub(sub_00, value_00, value_00, MPC_ROUNDING);
    mpc_t sub_01;
    mpc_init2(sub_01, MPC_PRECISION);
    mpc_sub(sub_01, value_01, value_01, MPC_ROUNDING);
    mpc_t sub_10;
    mpc_init2(sub_10, MPC_PRECISION);
    mpc_sub(sub_10, value_10, value_10, MPC_ROUNDING);
    mpc_t sub_11;
    mpc_init2(sub_11, MPC_PRECISION);
    mpc_sub(sub_11, value_11, value_11, MPC_ROUNDING);

    if(false) {
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)sub_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)sub_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)sub_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)sub_11, MPC_ROUNDING);
        putchar('\n');
    }

    test_assert( mtbdd_isleaf(dd_minus) );
    test_assert( mpc_compare( mtbdd_getvalue( dd_minus), (uint64_t)sub_00));
    test_assert( mpc_compare( mtbdd_getvalue( dd_minus), (uint64_t)sub_01));
    test_assert( mpc_compare( mtbdd_getvalue( dd_minus), (uint64_t)sub_10));
    test_assert( mpc_compare( mtbdd_getvalue( dd_minus), (uint64_t)sub_11));

    mpc_clear(sub_00);
    mpc_clear(sub_01);
    mpc_clear(sub_10);
    mpc_clear(sub_11);

    mpc_clear(value_00);
    mpc_clear(value_01);
    mpc_clear(value_10);
    mpc_clear(value_11);

    return 0;
}

int
test_mtbdd_arithmic_min_max_functions_complex()
{
    //
    // MTBDD dd1 =
    //
    //                         x1
    //
    //            x2                       x2
    //
    //  0.25+1.25i  0.75+1.75i    0.35+1.35i  0.65+1.65i
    //
    // The magnitudes are circa 1.272 1.905 1.395 1.774
    // 
    // MTBDD dd2 = 
    //
    //             x1
    //
    //      1.25 + 0.75i
    //
    // The magnitude is 1.457
    //
    // So, min(dd1,dd2) should result in
    // 
    //                     x1
    //            x2               x2
    //
    //  0.25+1.25i    1.25+0.75i    0.35+1.35i
    //

    MTBDD dd1, dd2;

    // Set the terminal leafs of dd1
    mpc_t value_00;
    mpc_assign(value_00, 0.25, 1.25);
    mpc_t value_01;
    mpc_assign(value_01, 0.75, 1.75);
    mpc_t value_10;
    mpc_assign(value_10, 0.35, 1.35);
    mpc_t value_11;
    mpc_assign(value_11, 0.65, 1.65);
    
    MTBDD index_leaf_00 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_00);
    MTBDD index_leaf_01 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_01);
    MTBDD index_leaf_10 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_10);
    MTBDD index_leaf_11 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_11);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x1_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x1_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    dd1 = index_root_node;

    // Set the terminal leaf of dd2
    mpc_t value_0;
    mpc_assign(value_0, 1.25, 0.75);
    
    MTBDD index_leaf_0 = mtbdd_makeleaf(MPC_TYPE, (size_t)value_0);

    // Make non-terminal nodes - middle layer, so variable x2
    index_x2 = 2;
    index_x1_low  = mtbdd_makenode(index_x2, index_leaf_0, index_leaf_0);
    index_x1_high = mtbdd_makenode(index_x2, index_leaf_0, index_leaf_0);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_root_node = mtbdd_makenode(index_x1, index_x1_low, index_x1_high);

    dd2 = index_root_node;

    // Compute minimum of a and b
    MTBDD dd_min = mpc_min(dd1, dd2);
    MTBDD dd_max = mpc_max(dd1, dd2);

    // Compute the answers
    mpc_t min_00;
    mpc_init2(min_00, MPC_PRECISION);
    mpc_minimum_abs(min_00, value_00, value_0);
    mpc_t min_01;
    mpc_init2(min_01, MPC_PRECISION);
    mpc_minimum_abs(min_01, value_01, value_0);
    mpc_t min_10;
    mpc_init2(min_10, MPC_PRECISION);
    mpc_minimum_abs(min_10, value_10, value_0);
    mpc_t min_11;
    mpc_init2(min_11, MPC_PRECISION);
    mpc_minimum_abs(min_11, value_11, value_0);

    // Compute the answers
    mpc_t max_00;
    mpc_init2(max_00, MPC_PRECISION);
    mpc_maximum_abs(max_00, value_00, value_0);
    mpc_t max_01;
    mpc_init2(max_01, MPC_PRECISION);
    mpc_maximum_abs(max_01, value_01, value_0);
    mpc_t max_10;
    mpc_init2(max_10, MPC_PRECISION);
    mpc_maximum_abs(max_10, value_10, value_0);
    mpc_t max_11;
    mpc_init2(max_11, MPC_PRECISION);
    mpc_maximum_abs(max_11, value_11, value_0);

    if(false) {
        printf("dd1 values \n");
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)value_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)value_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)value_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)value_11, MPC_ROUNDING);
        putchar('\n');

        printf("dd1 \n");
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_getlow(   mtbdd_getlow( dd1 ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_gethigh(  mtbdd_getlow( dd1 ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_getlow(  mtbdd_gethigh( dd1 ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_gethigh( mtbdd_gethigh( dd1 ))), MPC_ROUNDING);
        putchar('\n');

        printf("dd2 \n");
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue( dd2 ), MPC_ROUNDING);
        putchar('\n');

        printf("correct minimum values \n");
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)max_00, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)max_01, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)max_10, MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)max_11, MPC_ROUNDING);
        putchar('\n');

        printf("00 = %d\n", mtbdd_isleaf(mtbdd_getlow( mtbdd_getlow( dd_max))));
        printf("00 = %d\n", mtbdd_isleaf(mtbdd_gethigh(mtbdd_getlow( dd_max))));
        printf("00 = %d\n", mtbdd_isleaf(mtbdd_getlow( mtbdd_gethigh(dd_max))));
        printf("00 = %d\n", mtbdd_isleaf(mtbdd_gethigh(mtbdd_gethigh(dd_max))));

        printf("computed minimum values \n");
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_getlow(   mtbdd_getlow( dd_max ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_gethigh(  mtbdd_getlow( dd_max ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_getlow(  mtbdd_gethigh( dd_max ))), MPC_ROUNDING);
        putchar('\n');
        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(mtbdd_gethigh( mtbdd_gethigh( dd_max ))), MPC_ROUNDING);
        putchar('\n');
    }

    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow( mtbdd_getlow(  dd_min ))), (uint64_t)min_00));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_getlow(  dd_min ))), (uint64_t)min_01));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow( mtbdd_gethigh( dd_min ))), (uint64_t)min_10));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_gethigh( dd_min ))), (uint64_t)min_11));

    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow( mtbdd_getlow(  dd_max ))), (uint64_t)max_00));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_getlow(  dd_max ))), (uint64_t)max_01));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_getlow( mtbdd_gethigh( dd_max ))), (uint64_t)max_10));
    test_assert( mpc_compare( mtbdd_getvalue(mtbdd_gethigh(mtbdd_gethigh( dd_max ))), (uint64_t)max_11));

    mpc_clear(value_00);
    mpc_clear(value_01);
    mpc_clear(value_10);
    mpc_clear(value_11);

    mpc_clear(value_0);

    mpc_clear(min_00);
    mpc_clear(min_01);
    mpc_clear(min_10);
    mpc_clear(min_11);

    mpc_clear(max_00);
    mpc_clear(max_01);
    mpc_clear(max_10);
    mpc_clear(max_11);

    return 0;
}

// Testing mtbdd abstract arithmic functions double.

int
test_mtbdd_abstract_plus_function_1_double()
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

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

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

    test_assert(mtbdd_set_count(var_set) == length_var_set);
    test_assert(var[0] == var_[0]);
    if (length_var_set > 1) test_assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    //FILE *out = fopen("..//Testing//Temporary//output_dd_plus_1.dot", "w");
    //mtbdd_fprintdot(out, dd_plus);
    //fclose(out);

    test_assert(mtbdd_getdouble(mtbdd_getlow(dd_plus))  == 0.6);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.5);

    return 0;
}

int
test_mtbdd_abstract_plus_function_2_double()
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

    if(false) {
        printf("index_leaf_00 = %"PRIu64" \n", index_leaf_00);
        printf("index_leaf_01 = %"PRIu64" \n", index_leaf_01);
        printf("index_leaf_10 = %"PRIu64" \n", index_leaf_10);
        printf("index_leaf_11 = %"PRIu64" \n", index_leaf_11);
    }

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    if(false) {
        printf("index_x0_low  = %"PRIu64" \n", index_x0_low);
        printf("index_x0_high = %"PRIu64" \n", index_x0_high);
    }

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    //printf("index_x0 = %"PRIu64" \n", index_x0);

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

    test_assert(mtbdd_set_count(var_set) == length_var_set);
    test_assert(var[0] == var_[0]);
    if (length_var_set > 1) test_assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    //FILE *out = fopen("..//Testing//Temporary//output_dd_plus_2.dot", "w");
    //mtbdd_fprintdot(out, dd_plus);
    //fclose(out);

    // Print all kinds of gets
    if(false) {
        printf("dd_plus       = %"PRIu64"\n", dd_plus);
        printf("getnumer      = %d \n", mtbdd_getnumer(dd_plus));
        printf("getdouble     = %lf\n", mtbdd_getdouble(dd_plus));
        printf("getvalue      = %"PRIu64"\n", mtbdd_getvalue(dd_plus));
        printf("getlow        = %"PRIu64"\n", mtbdd_getlow(dd_plus));
        printf("gethigh       = %"PRIu64"\n", mtbdd_gethigh(dd_plus));
        printf("getvar        = %d \n", mtbdd_getvar(dd_plus));  // index_x1 (index_x0)

        printf("getlow(low)   00 = %lf\n", mtbdd_getdouble( mtbdd_getlow(dd_plus)));
        printf("getlow(high)  10 = %lf\n", mtbdd_getdouble( mtbdd_gethigh(dd_plus)));
    }

    test_assert(mtbdd_getdouble(mtbdd_getlow(dd_plus))  == 1.0);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(dd_plus)) == 1.1);

    return 0;
}

int
test_mtbdd_abstract_plus_min_max_times_function_3_double()
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

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

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

    test_assert(mtbdd_set_count(var_set) == length_var_set);
    test_assert(var[0] == var_[0]);
    if (length_var_set > 1) test_assert(var[1] == var_[1]);

    // Compute abstract_plus(dd, var_set)
    MTBDD dd_plus = mtbdd_abstract_plus(dd, var_set);

    // Print dd_plus
    //FILE *out = fopen("..//Testing//Temporary//output_dd_plus_3.dot", "w");
    //mtbdd_fprintdot(out, dd_plus);
    //fclose(out);

    test_assert(mtbdd_getdouble(dd_plus)  == 2.1);

    // Other operations
    MTBDD dd_min = mtbdd_abstract_min(dd, var_set);
    test_assert(mtbdd_getdouble(dd_min) == 0.25);

    MTBDD dd_max = mtbdd_abstract_max(dd, var_set);
    test_assert(mtbdd_getdouble(dd_max) == 0.75);

    MTBDD dd_times = mtbdd_abstract_times(dd, var_set);
    test_assert(mtbdd_getdouble(dd_times) == 0.25 * 0.75 * 0.35 * 0.75 );

    return 0;
}

// Testing mtbdd and abstract arithmic functions double.

/**
 *  Test with a matrix M[row][col] times vector v[row] multiplication
 *  
 *  M v = (m00 m01) (v0) = (m00.v0 + m01.v1) = w
 *        (m10 m11) (v1)   (m10.v0 + m11.v1)
 *
 *  Place m10 on leaf_01 and m01 on leaf_10 for M, and m01 on leaf_01 / m01 on leaf_01 for M_
 *
 *  Take MTBDD M =
 *
 *                           x0
 *                      0          1
 *                           x1
 *    col               0          1
 *                x2                      x2
 *    row    0          1          0            1
 *        
 *       m00 = 0.25  m10 = 0.35   m01 = 0.75   m11 = 0.65
 *
 *    M(x1,x2) = m00 . |x1 . |x2 + m10 . |x1 . x2 + m01 . x1 . |x2 + m11 . x1 . x2
 *
 *
 *  Take MTBDD v = 
 *
 *                           x0
 *                      0          1
 *                           x1
 *    row               0          1
 *
 *                   v0 = 3.0   v1 = 2.0
 *
 *    v(x1) = v0 . |x1 + v1 . x1
 *
 *  AND operation: 
 *
 *  w(x1,x2) = M(x1,x2) . v(x1)
 *           = m00 . v0 . |x1 . |x2 + m00 . v1 . |x1 . x1 . |x2 + m10 . v0 . |x1 . x2 + m10 . v1 . |x1 . x1 . x2 + ...
 *           = m00 . v0 . |x1 . |x2 + m10 . v0 . |x1 . x2 + m01 . v1 . x1 . |x2 + m11 . v1 . x1 . x2
 *
 *  w(x2)    = (m00 . v0 + m01 . v1) . |x2 + (m10 . v0 + m11 . v1) . x2, after elimination of x1
 *
 *  So, M v = w(x2), matrix multiplication of 2 x 2 . 2
 *
 *  Algorithm:
 *
 *    w(x2) = mtbdd_and_abstract_plus( M(x1,x2), v(x1), var_set = {1} )
 *
 *  Expected result:
 *
 *                   x0
 *                0      1
 *                   x2
 *    row         0      1
 *
 *               w0     w1
 *
 *
 *    w0 = 0.25 x 3.0 + 0.75 x 2.0
 *    w1 = 0.35 x 3.0 + 0.65 x 2.0
 */
int
test_mtbdd_and_abstract_plus_function_double()
{
    //// Create decision diagram for M(x1,x2)
    MTBDD M;

    // Set the terminal leafs
    MTBDD index_leaf_00 = mtbdd_double(0.25);
    MTBDD index_leaf_01 = mtbdd_double(0.35);
    MTBDD index_leaf_10 = mtbdd_double(0.75);
    MTBDD index_leaf_11 = mtbdd_double(0.65);

    // Make non-terminal nodes - middle layer, so variable x2
    uint32_t index_x2 = 2;
    MTBDD index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    MTBDD index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    uint32_t index_x1 = 1;
    MTBDD index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    M = index_x0;

    // Create decision diagram for M_(x1,x2)
    MTBDD M_;

    // Set the terminal leafs
    index_leaf_00 = mtbdd_double(0.25);
    index_leaf_01 = mtbdd_double(0.75);
    index_leaf_10 = mtbdd_double(0.35);
    index_leaf_11 = mtbdd_double(0.65);

    // Make non-terminal nodes - middle layer, so variable x2
    index_x2 = 2;
    index_x0_low  = mtbdd_makenode(index_x2, index_leaf_00, index_leaf_01);
    index_x0_high = mtbdd_makenode(index_x2, index_leaf_10, index_leaf_11);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_x0 = mtbdd_makenode(index_x1, index_x0_low, index_x0_high);

    M_ = index_x0;


    //// Create decision diagram v(x1)
    MTBDD v;

    // Set the terminal leafs
    MTBDD index_leaf_0 = mtbdd_double(3.0);
    MTBDD index_leaf_1 = mtbdd_double(2.0);

    // Make root node (= non terminal node) - top layer, so variable x1
    index_x1 = 1;
    index_x0 = mtbdd_makenode(index_x1, index_leaf_0, index_leaf_1);

    v = index_x0;


    //// Calculate M v = w

    // Prepare variable set to be removed from the v
    size_t length_var_set = 1;
    uint32_t var[length_var_set];
    var[0] = 1;
    
    // Test the mtbdd var_set to array and reverse function
    MTBDD var_set = mtbdd_set_from_array(var, length_var_set);

    // Compute and_abstract_plus(dd1, dd2, var_set)
    MTBDD w = mtbdd_and_abstract_plus(M, v, var_set);

    // Print w
    //FILE *out = fopen("..//Testing//Temporary//output_and_abstract_plus_w.dot", "w");
    //mtbdd_fprintdot(out, w);
    //fclose(out);

    double w0 = mtbdd_getdouble(mtbdd_getlow(w));
    double w1 = mtbdd_getdouble(mtbdd_gethigh(w));

    test_assert(mtbdd_getvar(w) == 2);
    test_assert(w0 == 0.25 * 3.0 + 0.75 * 2.0);
    test_assert(w1 == 0.35 * 3.0 + 0.65 * 2.0);

    //// Calculate M M v = M w = w_, the w has a node in x2 not x1, so use M_ because m10 <-> m01 

    // Prepare variable set to be removed from the v
    var[0] = 2;
    
    // Perpare the var_set
    var_set = mtbdd_set_from_array(var, length_var_set);
    MTBDD w_ = mtbdd_and_abstract_plus(M_, w, var_set);

    // Print w
    //out = fopen("..//Testing//Temporary//output_and_abstract_plus_mw.dot", "w");
    //mtbdd_fprintdot(out, w_);
    //fclose(out);

    test_assert(mtbdd_getvar(w_) == 1);
    test_assert(mtbdd_getdouble(mtbdd_getlow(w_))  == 0.25 * w0 + 0.75 * w1);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(w_)) == 0.35 * w0 + 0.65 * w1);

    return 0;
}

// Testing vector and matrix array conversion functions double.

int
test_vector_array_to_mtbdd_double()
{
    //
    //  From v = (1.0  2.0) make a matrix  V = (1.0 2.0)
    //                                         (1.0 2.0)
    //
    // mtbdd_to_vector_array is implemented with mtbdd_to_matrix_array
    //
    // So, first this last function is tested for conversion of vector_arrays
    //

    int n = 1; // Vector has size 2^n

    MTBDD v;
    
    // Fill dd column wise oriented
    v = mtbdd_makenode(0, 
        mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(2.0)),
        mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(2.0)));

    MatArr_t **W_arr = NULL;
    allocate_matrix_array(&W_arr, n);
    mtbdd_to_matrix_array(v, n, ALTERNATE_ROW_FIRST_WISE_MODE, W_arr);
    print_matrix_array(W_arr, n);

    test_assert(W_arr[0][0] == 1.0);
    test_assert(W_arr[0][1] == 1.0);
    test_assert(W_arr[1][0] == 2.0);
    test_assert(W_arr[1][1] == 2.0);

    // Fill dd column wise oriented
    v = mtbdd_makenode(0, 
        mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(1.0)),
        mtbdd_makenode(1, mtbdd_double(2.0), mtbdd_double(2.0)));

    allocate_matrix_array(&W_arr, n);
    mtbdd_to_matrix_array(v, n, ALTERNATE_ROW_FIRST_WISE_MODE, W_arr);
    print_matrix_array(W_arr, n);

    test_assert(W_arr[0][0] == 1.0);
    test_assert(W_arr[0][1] == 1.0);
    test_assert(W_arr[1][0] == 2.0);
    test_assert(W_arr[1][1] == 2.0);


    VecArr_t v_arr[(1 << n)];
    v_arr[0] = 0.0; v_arr[1] = 0.0;

    mtbdd_to_vector_array(v, n, ALTERNATE_ROW_FIRST_WISE_MODE, v_arr);

    print_vector_array(v_arr, n);

    test_assert(v_arr[0] == 1.0);
    test_assert(v_arr[1] == 2.0);

    v = vector_array_to_mtbdd(v_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    mtbdd_to_vector_array(v, n, ALTERNATE_ROW_FIRST_WISE_MODE, v_arr);

    print_vector_array(v_arr, n);

    test_assert(v_arr[0] == 1.0);
    test_assert(v_arr[1] == 2.0);

    return 0;
}

int
test_matrix_array_to_mtbdd_double()
{
    //
    //  K (x) L = M
    //
    //  K = (1.0  3.0)   L = (1.0  0.5)   M = (1.0 x L  3.0 x L)
    //      (2.0  1.0)       (0.5  1.0)       (2.0 x L  1.0 x L)
    //
    //  M = (1.0 0.5 3.0 1.5)
    //      (0.5 1.0 1.5 3.0)
    //      (2.0 1.0 1.0 0.5)
    //      (1.0 2.0 0.5 1.0)
    //

    int n = 1; // Matrix K and L have 2^n x 2^n size

    MTBDD K;
    
    // Fill both dd's row wise oriented
    K = mtbdd_makenode(0, 
        mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(3.0)),
        mtbdd_makenode(1, mtbdd_double(2.0), mtbdd_double(1.0)));

    MatArr_t **K_arr = NULL;
    test_assert(allocate_matrix_array(&K_arr, n) == 0);

    mtbdd_to_matrix_array(K, n, ALTERNATE_ROW_FIRST_WISE_MODE, K_arr);

    print_matrix_array(K_arr, n);

    test_assert(K_arr[0][0] == 1.0);
    test_assert(K_arr[0][1] == 3.0);
    test_assert(K_arr[1][0] == 2.0);
    test_assert(K_arr[1][1] == 1.0);

    K = matrix_array_to_mtbdd(K_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    mtbdd_to_matrix_array(K, n, ALTERNATE_ROW_FIRST_WISE_MODE, K_arr);

    print_matrix_array(K_arr, n);

    test_assert(K_arr[0][0] == 1.0);
    test_assert(K_arr[0][1] == 3.0);
    test_assert(K_arr[1][0] == 2.0);
    test_assert(K_arr[1][1] == 1.0);

    free_matrix_array(K_arr, n);

    return 0;
}

/**
 * K (x) L = M
 *
 * K = (1.0  3.0)   L = (1.0  0.5)   M = (1.0 x L  3.0 x L)
 *     (2.0  1.0)       (0.5  1.0)       (2.0 x L  1.0 x L)
 *
 * M = (1.0 0.5 3.0 1.5)
 *     (0.5 1.0 1.5 3.0)
 *     (2.0 1.0 1.0 0.5)
 *     (1.0 2.0 0.5 1.0)
 *
 */
int
test_mtbdd_to_matrix_array_double()
{
    int n = 2;

    MTBDD K, L, M;
    
    // Fill both dd's row wise oriented: f(col0, row0) = (K00 K01 K10 K11) = K_arr[row0][col0]
    K = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(3.0)),
                          mtbdd_makenode(1, mtbdd_double(2.0), mtbdd_double(1.0)));

    MatArr_t **K_arr = NULL;
    test_assert(allocate_matrix_array(&K_arr, 1) == 0);

    mtbdd_to_matrix_array(K, 1, ALTERNATE_ROW_FIRST_WISE_MODE, K_arr);

    print_matrix_array(K_arr, 1);

    test_assert(K_arr[0][0] == 1.0);
    test_assert(K_arr[0][1] == 3.0);
    test_assert(K_arr[1][0] == 2.0);
    test_assert(K_arr[1][1] == 1.0);

    free_matrix_array(K_arr, 1);

    L = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(0.5)),
                          mtbdd_makenode(1, mtbdd_double(0.5), mtbdd_double(1.0)));

    M = mtbdd_tensor_prod(K, L, n);

    MatArr_t **W_arr = NULL;
    test_assert(allocate_matrix_array(&W_arr, n) == 0);

    mtbdd_to_matrix_array(M, n, ALTERNATE_ROW_FIRST_WISE_MODE, W_arr);

    print_matrix_array(W_arr, n);

    test_assert(W_arr[0][0] == 1.0); // First row
    test_assert(W_arr[0][1] == 0.5);
    test_assert(W_arr[0][2] == 3.0);
    test_assert(W_arr[0][3] == 1.5);

    test_assert(W_arr[1][0] == 0.5); // Second row
    test_assert(W_arr[1][1] == 1.0);
    test_assert(W_arr[1][2] == 1.5);
    test_assert(W_arr[1][3] == 3.0);
    
    test_assert(W_arr[2][0] == 2.0); // Third row
    test_assert(W_arr[2][1] == 1.0);
    test_assert(W_arr[2][2] == 1.0);
    test_assert(W_arr[2][3] == 0.5);
    
    test_assert(W_arr[3][0] == 1.0); // Fourth row
    test_assert(W_arr[3][1] == 2.0);
    test_assert(W_arr[3][2] == 0.5);
    test_assert(W_arr[3][3] == 1.0);

    mtbdd_to_matrix_array(M, n, COLUMN_WISE_MODE, W_arr);

    print_matrix_array(W_arr, n);

    test_assert(W_arr[0][0] == 1.0); // Right Down
    test_assert(W_arr[0][1] == 0.5);
    test_assert(W_arr[0][2] == 0.5);
    test_assert(W_arr[0][3] == 1.0);

    test_assert(W_arr[1][0] == 3.0); // Right Up
    test_assert(W_arr[1][1] == 1.5);
    test_assert(W_arr[1][2] == 1.5);
    test_assert(W_arr[1][3] == 3.0);
    
    test_assert(W_arr[2][0] == 2.0); // Left Down
    test_assert(W_arr[2][1] == 1.0);
    test_assert(W_arr[2][2] == 1.0);
    test_assert(W_arr[2][3] == 2.0);
    
    test_assert(W_arr[3][0] == 1.0); // Left Up
    test_assert(W_arr[3][1] == 0.5);
    test_assert(W_arr[3][2] == 0.5);
    test_assert(W_arr[3][3] == 1.0);

    free_matrix_array(W_arr, n);

    return 0;
}

// Testing mtbdd kronecker functions double / complex

/**
 *  K (x) L = M
 *
 *  K = (1.0  3.0)   L = (1.0  0.5)   M = (1.0 x L  3.0 x L)
 *      (2.0  1.0)       (0.5  1.0)       (2.0 x L  1.0 x L)
 *  M = (1.0 0.5 3.0 1.5)
 *      (0.5 1.0 1.5 3.0)
 *      (2.0 1.0 1.0 0.5)
 *      (1.0 2.0 0.5 1.0) 
 */
int
test_mtbdd_matrix_kronecker_multiplication_double()
{
    MTBDD K, L, M;
    // Fill both dd's column wise oriented
    K = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(3.0)),
                          mtbdd_makenode(1, mtbdd_double(2.0), mtbdd_double(1.0)));
    L = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(0.5)),
                          mtbdd_makenode(1, mtbdd_double(0.5), mtbdd_double(1.0)));
    
    M = mtbdd_tensor_prod(K, L, 2); 
    
    // Read out column wise
    MTBDD M0000 = mtbdd_getlow(mtbdd_getlow(mtbdd_getlow(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0000));
    test_assert(mtbdd_getdouble(M0000) == 1.0);
    MTBDD M0001 = mtbdd_gethigh(mtbdd_getlow(mtbdd_getlow(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0001));
    test_assert(mtbdd_getdouble(M0001) == 0.5);
    MTBDD M0010 = mtbdd_getlow(mtbdd_gethigh(mtbdd_getlow(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0010));
    test_assert(mtbdd_getdouble(M0010) == 0.5);
    MTBDD M0011 = mtbdd_gethigh(mtbdd_gethigh(mtbdd_getlow(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0011));
    test_assert(mtbdd_getdouble(M0011) == 1.0);

    MTBDD M0100 = mtbdd_getlow(mtbdd_getlow(mtbdd_gethigh(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0100));
    test_assert(mtbdd_getdouble(M0100) == 3.0);
    MTBDD M0101 = mtbdd_gethigh(mtbdd_getlow(mtbdd_gethigh(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0101));
    test_assert(mtbdd_getdouble(M0101) == 1.5);
    MTBDD M0110 = mtbdd_getlow(mtbdd_gethigh(mtbdd_gethigh(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0110));
    test_assert(mtbdd_getdouble(M0110) == 1.5);
    MTBDD M0111 = mtbdd_gethigh(mtbdd_gethigh(mtbdd_gethigh(mtbdd_getlow(M))));
    test_assert(mtbdd_isleaf(M0111));
    test_assert(mtbdd_getdouble(M0111) == 3.0);

    MTBDD M1000 = mtbdd_getlow(mtbdd_getlow(mtbdd_getlow(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1000));
    test_assert(mtbdd_getdouble(M1000) == 2.0);
    MTBDD M1001 = mtbdd_gethigh(mtbdd_getlow(mtbdd_getlow(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1001));
    test_assert(mtbdd_getdouble(M1001) == 1.0);
    MTBDD M1010 = mtbdd_getlow(mtbdd_gethigh(mtbdd_getlow(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1010));
    test_assert(mtbdd_getdouble(M1010) == 1.0);
    MTBDD M1011 = mtbdd_gethigh(mtbdd_gethigh(mtbdd_getlow(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1011));
    test_assert(mtbdd_getdouble(M1011) == 2.0);

    MTBDD M1100 = mtbdd_getlow(mtbdd_getlow(mtbdd_gethigh(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1100));
    test_assert(mtbdd_getdouble(M1100) == 1.0);
    MTBDD M1101 = mtbdd_gethigh(mtbdd_getlow(mtbdd_gethigh(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1101));
    test_assert(mtbdd_getdouble(M1101) == 0.5);
    MTBDD M1110 = mtbdd_getlow(mtbdd_gethigh(mtbdd_gethigh(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1110));
    test_assert(mtbdd_getdouble(M1110) == 0.5);
    MTBDD M1111 = mtbdd_gethigh(mtbdd_gethigh(mtbdd_gethigh(mtbdd_gethigh(M))));
    test_assert(mtbdd_isleaf(M1111));
    test_assert(mtbdd_getdouble(M1111) == 1.0);

    return 0;
}

/**
 * K (x) L = M
 *
 * K = (K00  K01)   L = (L00  L01)   M = (K00 x L  K01 x L)
 *     (K10  K11)       (L10  L11)       (K10 x L  K11 x L)
 *
 * M = (K00.L00 K00.L01 K01.L00 K01.L01)
 *     (K00.L10 K00.L11 K01.L10 K01.L11)
 *     (K10.L00 K10.L01 K11.L00 K11.L01)
 *     (K10.L10 K10.L11 K11.L10 K11.L11)
 */
int
test_mtbdd_matrix_kronecker_multiplication_complex()
{
    int n = 2;

    mpc_t K00, K01, K10, K11; // row column
    mpc_t L00, L01, L10, L11;

    mpc_init2(K00, MPC_PRECISION);
    mpc_assign(K00, 3.0, 0.5);
    mpc_init2(K01, MPC_PRECISION);
    mpc_assign(K01, -1.0, 0.5);
    mpc_init2(K10, MPC_PRECISION);
    mpc_assign(K10, 2.0, -0.5);
    mpc_init2(K11, MPC_PRECISION);
    mpc_assign(K11, -1.5, -0.5);

    mpc_init2(L00, MPC_PRECISION);
    mpc_assign(L00, 2.0, 0.5);
    mpc_init2(L01, MPC_PRECISION);
    mpc_assign(L01, 3.0, 0.0);
    mpc_init2(L10, MPC_PRECISION);
    mpc_assign(L10, -1.0, 1.5);
    mpc_init2(L11, MPC_PRECISION);
    mpc_assign(L11, 4.0, -0.5);

    MTBDD K, L, M;
    // Fill both dd's column wise oriented
    K = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_makeleaf(MPC_TYPE, (size_t)K00), mtbdd_makeleaf(MPC_TYPE, (size_t)K01)),
                          mtbdd_makenode(1, mtbdd_makeleaf(MPC_TYPE, (size_t)K10), mtbdd_makeleaf(MPC_TYPE, (size_t)K11)));
    L = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_makeleaf(MPC_TYPE, (size_t)L00), mtbdd_makeleaf(MPC_TYPE, (size_t)L01)),
                          mtbdd_makenode(1, mtbdd_makeleaf(MPC_TYPE, (size_t)L10), mtbdd_makeleaf(MPC_TYPE, (size_t)L11)));
    
    M = mtbdd_tensor_prod(K, L, n); 

//
    mpc_ptr **K_arr = NULL;
    test_assert(allocate_matrix_array_mpc(&K_arr, 1) == 0);

    mtbdd_to_matrix_array_mpc(K, 1, ALTERNATE_ROW_FIRST_WISE_MODE, K_arr);

    print_matrix_array_mpc(K_arr, 1);

    test_assert( mpc_compare( (uint64_t)K_arr[0][0], (uint64_t)K00) );
    test_assert( mpc_compare( (uint64_t)K_arr[0][1], (uint64_t)K01) );
    test_assert( mpc_compare( (uint64_t)K_arr[1][0], (uint64_t)K10) );
    test_assert( mpc_compare( (uint64_t)K_arr[1][1], (uint64_t)K11) );

    free_matrix_array_mpc(K_arr, 1);

//
    mpc_ptr **W_arr = NULL;
    test_assert(allocate_matrix_array_mpc(&W_arr, n) == 0);

    mtbdd_to_matrix_array_mpc(M, n, COLUMN_WISE_MODE, W_arr);

    print_matrix_array_mpc(W_arr, n);

    // Row 0
    mpc_t W00, W01, W02, W03;
    mpc_multiplication(W00, K00, L00);
    mpc_multiplication(W01, K00, L01);
    mpc_multiplication(W02, K00, L10);
    mpc_multiplication(W03, K00, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[0][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][3], (uint64_t)W03 ) );

    // Row 1
    mpc_multiplication(W00, K01, L00);
    mpc_multiplication(W01, K01, L01);
    mpc_multiplication(W02, K01, L10);
    mpc_multiplication(W03, K01, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[1][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][3], (uint64_t)W03 ) );

    // Row 2
    mpc_multiplication(W00, K10, L00);
    mpc_multiplication(W01, K10, L01);
    mpc_multiplication(W02, K10, L10);
    mpc_multiplication(W03, K10, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[2][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[2][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[2][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[2][3], (uint64_t)W03 ) );

    // Row 3
    mpc_multiplication(W00, K11, L00);
    mpc_multiplication(W01, K11, L01);
    mpc_multiplication(W02, K11, L10);
    mpc_multiplication(W03, K11, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[3][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][3], (uint64_t)W03 ) );

    mtbdd_to_matrix_array_mpc(M, n, ALTERNATE_ROW_FIRST_WISE_MODE, W_arr);

    print_matrix_array_mpc(W_arr, n);

    // Left Up
    mpc_multiplication(W00, K00, L00);
    mpc_multiplication(W01, K00, L01);
    mpc_multiplication(W02, K00, L10);
    mpc_multiplication(W03, K00, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[0][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][0], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][1], (uint64_t)W03 ) );

    // Right Up
    mpc_multiplication(W00, K01, L00);
    mpc_multiplication(W01, K01, L01);
    mpc_multiplication(W02, K01, L10);
    mpc_multiplication(W03, K01, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[0][2], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][3], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][3], (uint64_t)W03 ) );

    // Left Down
    mpc_multiplication(W00, K10, L00);
    mpc_multiplication(W01, K10, L01);
    mpc_multiplication(W02, K10, L10);
    mpc_multiplication(W03, K10, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[2][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[2][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][0], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][1], (uint64_t)W03 ) );

    // Right Down
    mpc_multiplication(W00, K11, L00);
    mpc_multiplication(W01, K11, L01);
    mpc_multiplication(W02, K11, L10);
    mpc_multiplication(W03, K11, L11);

    test_assert( mpc_compare( (uint64_t)W_arr[2][2], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[2][3], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][2], (uint64_t)W02 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[3][3], (uint64_t)W03 ) );

    free_matrix_array_mpc(W_arr, n);

    mpc_clear(K00);
    mpc_clear(K01);
    mpc_clear(K10);
    mpc_clear(K11);

    mpc_clear(L00);
    mpc_clear(L01);
    mpc_clear(L10);
    mpc_clear(L11);

    mpc_clear(W00);
    mpc_clear(W01);
    mpc_clear(W02);
    mpc_clear(W03);

    return 0;
}

// Testing supporting functions for vector matrix multiplication double / complex.

/** 
 *  Decrease var numbers from 3 -> 0
 * 
 *             x3
 *        x4        x4
 *     1.0  2.0  3.0  4.0
 * 
 * 
 *  Expected output:
 * 
 *             x0
 *        x1        x1
 *     1.0  2.0  3.0  4.0
 */
int
test_renumber_variables_double()
{
    MTBDD M = mtbdd_makenode(3, 
            mtbdd_makenode(4, mtbdd_double(1.0), mtbdd_double(2.0)),
            mtbdd_makenode(4, mtbdd_double(3.0), mtbdd_double(4.0))
        );

    // Decrease the var numbers
    MTBDD W = mtbdd_renumber_variables(M, 0);

    test_assert(mtbdd_getvar(W) == 0);
    test_assert(mtbdd_getvar(mtbdd_getlow(W)) == 1);

    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(W)))  == 1.0);
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(W))) == 3.0);

    // Increase the var numbers
    W = mtbdd_renumber_variables(M, 6);

    test_assert(mtbdd_getvar(W) == 6);
    test_assert(mtbdd_getvar(mtbdd_getlow(W)) == 7);

    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_getlow(W)))  == 1.0);
    test_assert(mtbdd_getdouble(mtbdd_getlow(mtbdd_gethigh(W))) == 3.0);

    return 0;
}

int
test_determine_top_var_and_leafcount_double()
{
    // Test with identical leafvalues

    MTBDD M = MTBDD_ZERO;

    MTBDD node1 = mtbdd_makenode(4, mtbdd_double(3.0), mtbdd_double(3.0));
    MTBDD node2 = mtbdd_makenode(4, mtbdd_double(3.0), mtbdd_double(3.0));
    M = mtbdd_makenode(3, node1, node2);

    int maxvar = -1;
    int minvar = 100;
    int leafcount = 0;
    determine_top_var_and_leafcount(M, &minvar, &maxvar, &leafcount);

    test_assert(maxvar == -1);
    test_assert(minvar == 100);
    test_assert(leafcount == 1);

    // Test with different leafvalues

    node1 = mtbdd_makenode(4, mtbdd_double(3.0), mtbdd_double(1.0));
    node2 = mtbdd_makenode(4, mtbdd_double(1.0), mtbdd_double(3.0));
    M = mtbdd_makenode(3, node1, node2);

    maxvar = -1;
    minvar = 100;
    leafcount = 0;
    determine_top_var_and_leafcount(M, &minvar, &maxvar, &leafcount);

    test_assert(maxvar == 4);
    test_assert(minvar == 3);
    test_assert(leafcount == 4); // Should be two!

    return 0;
}

/**
 *  Test with 
 * 
 * M =  
 *      |
 *      x1
 *      |
 *      x2 -
 *      |
 */
int
test_mtbdd_get_children_of_var_double()
{
    double value1 = 3.0;
    double value2 = 1.0;

    MTBDD node1 = mtbdd_makenode(2, mtbdd_double(value1), mtbdd_double(value2));
    MTBDD node2 = mtbdd_makenode(2, mtbdd_double(value2), mtbdd_double(value1));
    MTBDD M = mtbdd_makenode(1, node1, node2);

    MTBDD M_low = M;
    MTBDD M_high = M;
    mtbdd_get_children_of_var(M, &M_low, &M_high, 1);

    test_assert(mtbdd_getvar(M_low) == 2);
    test_assert(mtbdd_getvar(M_high) == 2);

    test_assert(mtbdd_getdouble(mtbdd_getlow(M_low)) == value1);
    test_assert(mtbdd_getdouble(mtbdd_gethigh(M_high)) == value1);

    return 0;
}

// Testing mtbdd matrix vector multiplication functions double / complex.

/**
 *  M . v = w, M: 2^n x 2^n, v: 2^n x 1 (n_row x n_col)
 */
int
test_mtbdd_matrix_vector_multiplication_double()
{
    int n = 1;

    MatArr_t **M_arr = NULL;
    allocate_matrix_array(&M_arr, n);

    // even variables index rows
    // odd variables index columns
    M_arr[0][0] = -2.2; M_arr[0][1] =  1.2;
    M_arr[1][0] =  2.4; M_arr[1][1] = -1.4;
    MTBDD M = matrix_array_to_mtbdd(M_arr, n, ROW_WISE_MODE);

    VecArr_t v_arr[(1 << n)];
    // this is a column vector, so should be indexed by even variables (rows)
    v_arr[0] = 1.1; v_arr[1] = -2.2;
    MTBDD v = vector_array_to_mtbdd(v_arr, n, ROW_WISE_MODE);

    int currentvar = 0;
    MTBDD product = mtbdd_matvec_mult(M, v, 2*n, currentvar);

    MatArr_t w_arr[2];
    // TODO: fix inconsistency between this COLUMN_WISE_MODE and ROW_WISE_MODE used for input
    mtbdd_to_vector_array(product, n, COLUMN_WISE_MODE, w_arr);

    test_assert(w_arr[0] == M_arr[0][0] * v_arr[0] + M_arr[0][1] * v_arr[1]);
    test_assert(w_arr[1] == M_arr[1][0] * v_arr[0] + M_arr[1][1] * v_arr[1]);

    free_matrix_array(M_arr, n);

    return 0;
}

/** 
*  M . v = w, M: 2^n x 2^n, v: 2^n x 1 (n_row x n_col)
*
* This alternative matrix vector multiplication function 
* keeps the vector index orientation intact (v and w column wise). 
*
* Only n < 2.
*/
int
test_mtbdd_matrix_vector_multiplication_alt_double()
{
    int n = 1;

    MatArr_t **M_arr = NULL;
    allocate_matrix_array(&M_arr, n);

    M_arr[0][0] = -2.2; M_arr[0][1] =  1.2;
    M_arr[1][0] =  2.4; M_arr[1][1] = -1.4;
    MTBDD M = matrix_array_to_mtbdd(M_arr, n, ROW_WISE_MODE);

    VecArr_t v_arr[(1 << n)];
    v_arr[0] = 1.1; v_arr[1] = -2.2;
    MTBDD v = vector_array_to_mtbdd(v_arr, n, COLUMN_WISE_MODE);

    MTBDD product = mtbdd_matvec_mult_alt(M, v, n);

    MatArr_t w_arr[2];
    mtbdd_to_vector_array(product, n, COLUMN_WISE_MODE, w_arr);

    test_assert(w_arr[0] == M_arr[0][0] * v_arr[0] + M_arr[0][1] * v_arr[1]);
    test_assert(w_arr[1] == M_arr[1][0] * v_arr[0] + M_arr[1][1] * v_arr[1]);

    free_matrix_array(M_arr, n);

    return 0;
}

/**
 *  K . L = W, K: 2^n x 2^n, L, W: 2^n x 1 (n_row x n_col)
 */
int
test_mtbdd_matrix_vector_multiplication_complex()
{
    int n = 1;

    mpc_ptr **K_arr = NULL;
    allocate_matrix_array_mpc(&K_arr, n);

    mpc_t K00, K01, K10, K11; // row column
    mpc_t L00, L10;

    mpc_init2(K00, MPC_PRECISION);
    mpc_assign(K00, 3.0, 0.5);
    mpc_init2(K01, MPC_PRECISION);
    mpc_assign(K01, -1.0, 0.5);
    mpc_init2(K10, MPC_PRECISION);
    mpc_assign(K10, 2.0, -0.5);
    mpc_init2(K11, MPC_PRECISION);
    mpc_assign(K11, -1.5, -0.5);

    mpc_init2(L00, MPC_PRECISION);
    mpc_assign(L00, 2.0, 0.5);
    mpc_init2(L10, MPC_PRECISION);
    mpc_assign(L10, 3.0, 0.0);

    // even variables index rows
    // odd variables index columns
    K_arr[0][0] = K00; K_arr[0][1] = K01;
    K_arr[1][0] = K10; K_arr[1][1] = K11;
    MTBDD K = matrix_array_to_mtbdd_mpc(K_arr, n, ROW_WISE_MODE);

    mpc_ptr v_arr[(1 << n)];
    // this is a column vector, so should be indexed by even variables (rows)
    v_arr[0] = L00; v_arr[1] = L10;
    MTBDD v = vector_array_to_mtbdd_mpc(v_arr, n, ROW_WISE_MODE);

    int currentvar = 0;
    MTBDD product = mtbdd_matvec_mult(K, v, 2*n, currentvar);

    mpc_ptr W_arr[2];
    // TODO: fix inconsistency between this COLUMN_WISE_MODE and ROW_WISE_MODE used for input
    mtbdd_to_vector_array_mpc(product, n, COLUMN_WISE_MODE, W_arr);

    mpc_t X00, X10; // X[row, col]
    mpc_t Y00, Y10;
    mpc_t W00, W10;

    mpc_multiplication(X00, K00, L00);
    mpc_multiplication(X10, K10, L00);

    mpc_multiplication(Y00, K01, L10);
    mpc_multiplication(Y10, K11, L10);

    mpc_addition(W00, X00, Y00);
    mpc_addition(W10, X10, Y10);

    test_assert( mpc_compare( (uint64_t)W_arr[0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1], (uint64_t)W10 ) );

    free_matrix_array_mpc(K_arr, n);

    mpc_clear(K00);
    mpc_clear(K01);
    mpc_clear(K10);
    mpc_clear(K11);

    mpc_clear(L00);
    mpc_clear(L10);

    mpc_clear(X00);
    mpc_clear(X10);

    mpc_clear(Y00);
    mpc_clear(Y10);

    mpc_clear(W00);
    mpc_clear(W10);

    return 0;
}

// Testing mtbdd matrix matrix multiplication functions double / complex.

/** 
 *  K . M = W, M: 2^n x 2^n, L: 2^n x 2^n, W: 2^n x 2^n
 */ 
int
test_mtbdd_matrix_matrix_multiplication_1_double()
{
    int n = 1;

    MatArr_t **K_arr = NULL;
    allocate_matrix_array(&K_arr, n);

    K_arr[0][0] =  1.0; K_arr[0][1] = 0.0; 
    K_arr[1][0] =  1.0; K_arr[1][1] = 1.0;
    MTBDD K = matrix_array_to_mtbdd(K_arr, n, ROW_WISE_MODE);

    MatArr_t **M_arr = NULL;
    allocate_matrix_array(&M_arr, n);

    M_arr[0][0] =  1.0; M_arr[0][1] =  0.0;
    M_arr[1][0] =  0.0; M_arr[1][1] =  2.0;
    MTBDD M = matrix_array_to_mtbdd(M_arr, n, ROW_WISE_MODE);

    int currentvar = 0;
    MTBDD product = mtbdd_matmat_mult(K, M, 2*n, currentvar);

    MatArr_t **W_arr = NULL;
    allocate_matrix_array(&W_arr, n);
    mtbdd_to_matrix_array(product, n, COLUMN_WISE_MODE, W_arr);

    print_matrix_array(K_arr, n);
    print_matrix_array(M_arr, n);
    print_matrix_array(W_arr, n);

    if(false) {
        printf("W[0][0]: %lf %lf\n", W_arr[0][0], K_arr[0][0] * M_arr[0][0] + K_arr[0][1] * M_arr[1][0]);
        printf("W[0][1]: %lf %lf\n", W_arr[0][1], K_arr[0][0] * M_arr[0][1] + K_arr[0][1] * M_arr[1][1]);
        printf("W[1][0]: %lf %lf\n", W_arr[1][0], K_arr[1][0] * M_arr[0][0] + K_arr[1][1] * M_arr[1][0]);
        printf("W[1][1]: %lf %lf\n", W_arr[1][1], K_arr[1][0] * M_arr[0][1] + K_arr[1][1] * M_arr[1][1]);
    }

    test_assert(W_arr[0][0] == K_arr[0][0] * M_arr[0][0] + K_arr[0][1] * M_arr[1][0]);
    test_assert(W_arr[0][1] == K_arr[0][0] * M_arr[0][1] + K_arr[0][1] * M_arr[1][1]);

    test_assert(W_arr[1][0] == K_arr[1][0] * M_arr[0][0] + K_arr[1][1] * M_arr[1][0]);
    test_assert(W_arr[1][1] == K_arr[1][0] * M_arr[0][1] + K_arr[1][1] * M_arr[1][1]);

    free_matrix_array(K_arr, n);
    free_matrix_array(M_arr, n);
    free_matrix_array(W_arr, n);

    return 0;
}

/**
 *  K . M = W, M: 2^n x 2^n, L: 2^n x 2^n, W: 2^n x 2^n
 */
int
test_mtbdd_matrix_matrix_multiplication_2_double()
{
    int n = 1;

    MatArr_t **K_arr = NULL;
    allocate_matrix_array(&K_arr, n);

    K_arr[0][0] =  1.0; K_arr[0][1] = 2.0; 
    K_arr[1][0] =  1.0; K_arr[1][1] = 3.0;
    MTBDD K = matrix_array_to_mtbdd(K_arr, n, ROW_WISE_MODE);

    MatArr_t **M_arr = NULL;
    allocate_matrix_array(&M_arr, n);

    M_arr[0][0] =  1.0; M_arr[0][1] =  2.0;
    M_arr[1][0] =  1.0; M_arr[1][1] =  2.0;
    MTBDD M = matrix_array_to_mtbdd(M_arr, n, ROW_WISE_MODE);

    int currentvar = 0;
    MTBDD product = mtbdd_matmat_mult(K, M, 2*n, currentvar);

    MatArr_t **W_arr = NULL;
    allocate_matrix_array(&W_arr, n);
    mtbdd_to_matrix_array(product, n, COLUMN_WISE_MODE, W_arr);

    print_matrix_array(K_arr, n);
    print_matrix_array(M_arr, n);
    print_matrix_array(W_arr, n);

    if(false) {
        printf("W[0][0]: %lf %lf\n", W_arr[0][0], K_arr[0][0] * M_arr[0][0] + K_arr[0][1] * M_arr[1][0]);
        printf("W[0][1]: %lf %lf\n", W_arr[0][1], K_arr[0][0] * M_arr[0][1] + K_arr[0][1] * M_arr[1][1]);
        printf("W[1][0]: %lf %lf\n", W_arr[1][0], K_arr[1][0] * M_arr[0][0] + K_arr[1][1] * M_arr[1][0]);
        printf("W[1][1]: %lf %lf\n", W_arr[1][1], K_arr[1][0] * M_arr[0][1] + K_arr[1][1] * M_arr[1][1]);
    }

    test_assert(W_arr[0][0] == K_arr[0][0] * M_arr[0][0] + K_arr[0][1] * M_arr[1][0]);
    test_assert(W_arr[0][1] == K_arr[0][0] * M_arr[0][1] + K_arr[0][1] * M_arr[1][1]);

    test_assert(W_arr[1][0] == K_arr[1][0] * M_arr[0][0] + K_arr[1][1] * M_arr[1][0]);
    test_assert(W_arr[1][1] == K_arr[1][0] * M_arr[0][1] + K_arr[1][1] * M_arr[1][1]);

    free_matrix_array(K_arr, n);
    free_matrix_array(M_arr, n);
    free_matrix_array(W_arr, n);

    return 0;
}

/**
 *  K (x) L = M
 * 
 *  K = (1.0  3.0)   L = (1.0  0.5)   M = (1.0 x L  3.0 x L)
 *      (2.0  1.0)       (0.5  1.0)       (2.0 x L  1.0 x L)
 * 
 *  M = (1.0 0.5 3.0 1.5)
 *      (0.5 1.0 1.5 3.0)
 *      (2.0 1.0 1.0 0.5)
 *      (1.0 2.0 0.5 1.0)
 * 
 *  W = M.M
 * 
 *  W = ()
 * 
 *  Test is to check if the multiplication is correct for 4 x 4 sized matrices.
 */
int
test_matrix_matrix_multiplication_4x4_double()
{
    // Compose matrix M with expansive Kronecker product
    int n = 2;

    MTBDD K, L, M1, M2, W;
    
    // Fill both dd's column wise oriented
    K = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(3.0)),
                          mtbdd_makenode(1, mtbdd_double(2.0), mtbdd_double(4.0)));

    L = mtbdd_makenode(0, mtbdd_makenode(1, mtbdd_double(1.0), mtbdd_double(2.0)),
                          mtbdd_makenode(1, mtbdd_double(0.5), mtbdd_double(1.5)));

    M1 = mtbdd_tensor_prod(K, L, n);
    M2 = mtbdd_tensor_prod(K, L, n);

    test_assert(M1 == M2);

    // Calculate W = M1.M2 = M.M
    int currentvar = 0;
    W = mtbdd_matmat_mult(M1, M2, 2*n, currentvar);

    // Evaluate with matrix arrays

    MatArr_t **W_ = NULL;
    allocate_matrix_array(&W_, n);
    mtbdd_to_matrix_array(W, n, ALTERNATE_ROW_FIRST_WISE_MODE, W_);

    MatArr_t **M_ = NULL;
    allocate_matrix_array(&M_, n);
    mtbdd_to_matrix_array(M1, n, ALTERNATE_ROW_FIRST_WISE_MODE, M_);

    print_matrix_array(W_, n);
    print_matrix_array(M_, n);

    for(int col=0; col<(1<<n); col++) { // column in second matrix
        for(int row=0; row<(1<<n); row++) { // row in first matrix
            double element = 0.0;
            for(int iterator=0; iterator<(1<<n); iterator++)
                element += M_[row][iterator] * M_[iterator][col];
            test_assert(W_[row][col] == element);
        }
    }

    free_matrix_array(W_, n);
    free_matrix_array(M_, n);

    return 0;
}

/** 
 *  K . M = W, M: n x n, L: 2^n x 2^n, W: 2^n x 2^n
 */ 
int
test_mtbdd_matrix_matrix_multiplication_alt_double()
{
    int n = 1;

    MatArr_t **K_arr = NULL;
    allocate_matrix_array(&K_arr, n);

    K_arr[0][0] =  1.0; K_arr[0][1] = 2.0; 
    K_arr[1][0] = -1.0; K_arr[1][1] = 3.0;
    MTBDD K = matrix_array_to_mtbdd(K_arr, n, ROW_WISE_MODE);

    MatArr_t **M_arr = NULL;
    allocate_matrix_array(&M_arr, n);

    M_arr[0][0] = -1.0; M_arr[0][1] =  2.0;
    M_arr[1][0] =  1.0; M_arr[1][1] = -2.0;
    MTBDD M = matrix_array_to_mtbdd(M_arr, n, COLUMN_WISE_MODE);

    MTBDD product = mtbdd_matmat_mult_alt(K, M, n);

    MatArr_t **W_arr = NULL;
    allocate_matrix_array(&W_arr, n);
    mtbdd_to_matrix_array(product, n, ROW_WISE_MODE, W_arr);

    print_matrix_array(K_arr, n);
    print_matrix_array(M_arr, n);
    print_matrix_array(W_arr, n);

    test_assert(W_arr[0][0] == K_arr[0][0] * M_arr[0][0] + K_arr[0][1] * M_arr[1][0]);
    test_assert(W_arr[0][1] == K_arr[0][0] * M_arr[0][1] + K_arr[0][1] * M_arr[1][1]);

    test_assert(W_arr[1][0] == K_arr[1][0] * M_arr[0][0] + K_arr[1][1] * M_arr[1][0]);
    test_assert(W_arr[1][1] == K_arr[1][0] * M_arr[0][1] + K_arr[1][1] * M_arr[1][1]);

    free_matrix_array(K_arr, n);
    free_matrix_array(M_arr, n);
    free_matrix_array(W_arr, n);

    return 0;
}

/**
 *  K . L = W, K: 2^n x 2^n, L: 2^n x 2^n, W: 2^n x 2^n
 */ 
int
test_mtbdd_matrix_matrix_multiplication_1_complex()
{
    int n = 1;

    mpc_ptr **K_arr = NULL;
    allocate_matrix_array_mpc(&K_arr, n);

    mpc_t K00, K01, K10, K11; // row column
    mpc_t L00, L01, L10, L11;

    mpc_init2(K00, MPC_PRECISION);
    mpc_assign(K00, 3.0, 0.5);
    mpc_init2(K01, MPC_PRECISION);
    mpc_assign(K01, -1.0, 0.5);
    mpc_init2(K10, MPC_PRECISION);
    mpc_assign(K10, 2.0, -0.5);
    mpc_init2(K11, MPC_PRECISION);
    mpc_assign(K11, -1.5, -0.5);

    mpc_init2(L00, MPC_PRECISION);
    mpc_assign(L00, 2.0, 0.5);
    mpc_init2(L01, MPC_PRECISION);
    mpc_assign(L01, 3.0, 0.0);
    mpc_init2(L10, MPC_PRECISION);
    mpc_assign(L10, -1.0, 1.5);
    mpc_init2(L11, MPC_PRECISION);
    mpc_assign(L11, 4.0, -0.5);

    K_arr[0][0] = K00; K_arr[0][1] = K01; 
    K_arr[1][0] = K10; K_arr[1][1] = K11;
    MTBDD K = matrix_array_to_mtbdd_mpc(K_arr, n, ROW_WISE_MODE);

    mpc_ptr **L_arr = NULL;
    allocate_matrix_array_mpc(&L_arr, n);

    L_arr[0][0] = L00; L_arr[0][1] = L01;
    L_arr[1][0] = L10; L_arr[1][1] = L11;
    MTBDD L = matrix_array_to_mtbdd_mpc(L_arr, n, ROW_WISE_MODE);

    int currentvar = 0;
    MTBDD product = mtbdd_matmat_mult(K, L, 2*n, currentvar);

    mpc_ptr **W_arr = NULL;
    allocate_matrix_array_mpc(&W_arr, n);
    mtbdd_to_matrix_array_mpc(product, n, COLUMN_WISE_MODE, W_arr);

    print_matrix_array_mpc(K_arr, n);
    print_matrix_array_mpc(L_arr, n);
    print_matrix_array_mpc(W_arr, n);

    mpc_t X00, X01, X10, X11;
    mpc_t Y00, Y01, Y10, Y11;
    mpc_t W00, W01, W10, W11;

    mpc_multiplication(X00, K00, L00);
    mpc_multiplication(X01, K00, L01);
    mpc_multiplication(X10, K10, L00);
    mpc_multiplication(X11, K10, L01);

    mpc_multiplication(Y00, K01, L10);
    mpc_multiplication(Y01, K01, L11);
    mpc_multiplication(Y10, K11, L10);
    mpc_multiplication(Y11, K11, L11);

    mpc_addition(W00, X00, Y00);
    mpc_addition(W01, X01, Y01);
    mpc_addition(W10, X10, Y10);
    mpc_addition(W11, X11, Y11);

    test_assert( mpc_compare( (uint64_t)W_arr[0][0], (uint64_t)W00 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[0][1], (uint64_t)W01 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][0], (uint64_t)W10 ) );
    test_assert( mpc_compare( (uint64_t)W_arr[1][1], (uint64_t)W11 ) );

    free_matrix_array_mpc(K_arr, n);
    free_matrix_array_mpc(L_arr, n);
    free_matrix_array_mpc(W_arr, n);

    mpc_clear(K00);
    mpc_clear(K01);
    mpc_clear(K10);
    mpc_clear(K11);

    mpc_clear(L00);
    mpc_clear(L01);
    mpc_clear(L10);
    mpc_clear(L11);

    mpc_clear(X00);
    mpc_clear(X01);
    mpc_clear(X10);
    mpc_clear(X11);

    mpc_clear(Y00);
    mpc_clear(Y01);
    mpc_clear(Y10);
    mpc_clear(Y11);

    mpc_clear(W00);
    mpc_clear(W01);
    mpc_clear(W10);
    mpc_clear(W11);

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

TASK_0(int, runtests)
{
    // We are not testing garbage collection
    sylvan_gc_disable();

    // Test 1
    printf("\nTesting mtbdd makenode and ithvar.\n");
    if (test_mtbdd_makenode_ithvar()) return 1;

    // Test 2
    printf("\nTesting mtbdd makeleaf, makenode, leaf type boolean, integer, double, complex.\n");
    if (test_mtbdd_makenodes_and_leafs_boolean()) return 1;
    if (test_mtbdd_makenodes_and_leafs_integer()) return 1;
    if (test_mtbdd_makenodes_and_leafs_double()) return 1;
    if (test_mtbdd_makenodes_and_leafs_complex()) return 1;

    // Test 3.1
    printf("\nTesting mtbdd arithmic functions double.\n");
    if (test_mtbdd_arithmic_functions_double()) return 1;
    
    // Test 3.2
    printf("\nTesting mtbdd arithmic functions complex.\n");
    if (test_mtbdd_arithmic_plus_sub_times_functions_complex()) return 1;
    if (test_mtbdd_arithmic_min_max_functions_complex()) return 1;

    // Test 4.1
    printf("\nTesting mtbdd abstract arithmic functions double.\n");
    if (test_mtbdd_abstract_plus_function_1_double()) return 1;
    if (test_mtbdd_abstract_plus_function_2_double()) return 1;
    if (test_mtbdd_abstract_plus_min_max_times_function_3_double()) return 1;

    // Test 4.2
    printf("\nTesting mtbdd abstract arithmic functions complex.\n");
    // if (test_mtbdd_abstract_plus_function_1_complex()) return 1;
    // if (test_mtbdd_abstract_plus_function_2_complex()) return 1;
    // if (test_mtbdd_abstract_plus_min_max_times_function_3_complex()) return 1;

    // Test 5.1
    printf("\nTesting mtbdd and abstract arithmic functions double.\n");
    if (test_mtbdd_and_abstract_plus_function_double()) return 1;

    // Test 5.2
    printf("\nTesting mtbdd and abstract arithmic functions complex.\n");
    // if (test_mtbdd_and_abstract_plus_function_complex()) return 1;

    // Test 6.1
    printf("\nTesting vector and matrix array conversion functions double.\n");
    if (test_vector_array_to_mtbdd_double()) return 1;
    if (test_matrix_array_to_mtbdd_double()) return 1;
    if (test_mtbdd_to_matrix_array_double()) return 1;

    // Test 6.2
    printf("\nTesting vector and matrix array conversion functions double.\n");
    // if (test_vector_array_to_mtbdd_complex()) return 1;
    // if (test_matrix_array_to_mtbdd_complex()) return 1;
    // if (test_mtbdd_to_matrix_array_complex()) return 1;

    // Test 7.1
    printf("\nTesting mtbdd kronecker functions double.\n");
    if (test_mtbdd_matrix_kronecker_multiplication_double()) return 1;

    // Test 7.2
    printf("\nTesting mtbdd kronecker functions complex.\n");
    if (test_mtbdd_matrix_kronecker_multiplication_complex()) return 1;

    // Test 8.1
    printf("\nTesting supporting functions for vector matrix multiplication double.\n");
    if (test_renumber_variables_double()) return 1;
    if (test_determine_top_var_and_leafcount_double()) return 1;
    if (test_mtbdd_get_children_of_var_double()) return 1;

    // Test 8.2
    printf("\nTesting supporting functions for vector matrix multiplication complex.\n");
    // if (test_renumber_variables_complex()) return 1;
    // if (test_determine_top_var_and_leafcount_complex()) return 1;
    // if (test_mtbdd_get_children_of_var_complex()) return 1;

    // Test 9.1
    printf("\nTesting mtbdd matrix vector multiplication functions double.\n");
    if (test_mtbdd_matrix_vector_multiplication_double()) return 1;
    if (test_mtbdd_matrix_vector_multiplication_alt_double()) return 1;

    // Test 9.2
    printf("\nTesting mtbdd matrix vector multiplication functions complex.\n");
    if (test_mtbdd_matrix_vector_multiplication_complex()) return 1;
    // if (test_mtbdd_matrix_vector_multiplication_alt_complex()) return 1;

    // Test 10.1
    printf("\nTesting mtbdd matrix matrix multiplication functions double.\n");
    if (test_mtbdd_matrix_matrix_multiplication_1_double()) return 1;
    if (test_mtbdd_matrix_matrix_multiplication_2_double()) return 1;
    if (test_matrix_matrix_multiplication_4x4_double()) return 1;
    if (test_mtbdd_matrix_matrix_multiplication_alt_double()) return 1;

    // Test 10.2
    printf("\nTesting mtbdd matrix matrix multiplication functions complex.\n");
    if (test_mtbdd_matrix_matrix_multiplication_1_complex()) return 1;
    // if (test_mtbdd_matrix_matrix_multiplication_2_complex()) return 1;
    // if (test_matrix_matrix_multiplication_4x4_complex()) return 1;
    // if (test_mtbdd_matrix_matrix_multiplication_alt_complex()) return 1;

    return 0;
}

int main()
{
    // Standard Lace initialization with 1 worker
    lace_start(1, 0);

    // Simple Sylvan initialization, also initialize BDD, MTBDD and LDD support
    sylvan_set_sizes(1LL<<20, 1LL<<20, 1LL<<16, 1LL<<16);
    sylvan_init_package(); 
    
    sylvan_init_mtbdd();

    printf("Mtbdd initialization complete.\n\n");

    uint32_t mpc_type = mpc_init();
    printf("Mtbdd mpc type initialization complete, mpc_type = %d.\n\n", mpc_type);

    test_assert(mpc_type == MPC_TYPE);

    int result = RUN(runtests);

    sylvan_quit();
    lace_stop();

    return result;
}



