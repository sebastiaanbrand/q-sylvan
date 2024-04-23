/*
 * Copyright 2024 System Verification Lab, LIACS, Leiden University
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

#include <qsylvan_simulator_mtbdd.h>
#include <sylvan.h>
#include <sylvan_int.h>
#include <sylvan_mpc.h>

/**
 * Convert zero state vector x[] = {0.0, 0.0, ..., 0.0} to MTBDD 
 * 
 * This results in an MTBDD with mpc complex zero leafs, deepness n
 * 
 *              x0
 *          x1      0.0+i0.0
 *  ...
 *          x(n-1)
 *  0.0+i0.0      0.0+i0.0
 *
 */
MTBDD
mtbdd_create_all_zero_state_double(BDDVAR n)
{
    bool x[n];
    for (BDDVAR k=0; k<n; k++) x[k] = 0;
    return mtbdd_create_basis_state_double(n, x);
}

MTBDD
mtbdd_create_all_zero_state_mpc(BDDVAR n)
{
    bool x[n];
    for (BDDVAR k=0; k<n; k++) x[k] = 0;
    return mtbdd_create_basis_state_mpc(n, x);
}

/**
 * 
 * Convert one state column vector s = (0.0, 1.0, ..., 0.0, 0.0) 
 * to MTBDD
 * 
 * Example:
 * 
 *   |101> = |(0,0,0,0,1,0,0,0) = |4+1> = |5>
 * 
 * Equivalent with:
 * 
 *   x[n] = {false, false, ..., false}
 * 
 * This results in an MTBDD with mpc complex leafs, deepness n,
 * so number of leafs 2 ** n.
 * 
 * The vars are even {0,2,4,...} to align linear algebra operations.
 * 
 *                          x0
 * 
 *             x2                        x2
 *  
 *        ....                         ....
 * 
 *   x(2**(n-1))  x(2**(n-1))  x(2**(n-1))  x(2**(n-1))
 * 
 *   0.0   0.0    0.0   0.0    1.0   0.0    0.0   0.0
 * 
 * This will be reduced to
 * 
 *                    x0
 * 
 *             x2            0.0
 *  
 *          /        0.0
 * 
 *   x(2**(n-1))  
 * 
 *   0.0      0.0 
 * 
 * We will build the reduced mtbdd directly starting at the bottom where we place a 1 leaf.
 * 
 *      n       highest var
 *      0       0
 *      1       n-1 = 0
 *      2       n+0 = 2
 *      3       n+1 = 4
 *      4       n+2 = 6 
 *         ...
 *      n       n+n-2 = 2n-2, n > 0
 * 
 */
MTBDD
mtbdd_create_basis_state_double(BDDVAR n, bool* x)
{
    if(n==0)
        return MTBDD_ZERO;

    uint32_t var = 2*n-2;

    double double_zero = 0.0;
    double double_one = 1.0;

    MTBDD zero = mtbdd_double(double_zero);
    MTBDD one = mtbdd_double(double_one);

    MTBDD node = one;

    //
    // Start with the least significant qubit first, x[n-1]
    //
    // Build a path from the bottom (place one leaf = 1.0) leaf to the root.
    //
    // If the qubit is 0 choose the low edge, otherwise the high edge
    //
    for(int i = (int)n - 1; i>=0; i--)
    {
        printf("x[%d] = %d, var = %d\n", i, x[i], var);

        if(x[i] == 0)
            node = mtbdd_makenode(var, node, zero);

        if(x[i] == 1)
            node = mtbdd_makenode(var, zero, node);

        // var of node always even
        var = var - 2;
    }

    return node;
}

MTBDD
mtbdd_create_basis_state_mpc(BDDVAR n, bool* x)
{
    if(n==0)
        return MTBDD_ZERO;

    uint32_t var = 2*n-2;

    mpc_t mpc_zero;
    mpc_init2(mpc_zero, MPC_PRECISION);
    mpc_assign(mpc_zero, 0.0, 0.0);

    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    MTBDD zero = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_zero);
    MTBDD one = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_one);

    MTBDD node = one;

    //
    // Start with the least significant qubit first, x[n-1]
    //
    // Build a path from the bottom (place one leaf = 1.0) leaf to the root.
    //
    // If the qubit is 0 choose the low edge, otherwise the high edge
    //
    for(int i = (int)n - 1; i>=0; i--)
    {
        printf("x[%d] = %d, var = %d\n", i, x[i], var);

        if(x[i] == 0)
            node = mtbdd_makenode(var, node, zero);

        if(x[i] == 1)
            node = mtbdd_makenode(var, zero, node);

        // var of node always even
        var = var - 2;
    }

    mpc_clear(mpc_zero);
    mpc_clear(mpc_one);

    return node;
}

/**
 * Create a single Qubit gate
 * 
 * 
 */
MTBDD
mtbdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid)
{
    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    MTBDD prev = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_one);
    for (int k = n-1; k >= 0; k--) {
        if ((unsigned int)k == t)
            prev = mtbdd_stack_matrix(k, gateid);
        else
            prev = mtbdd_stack_matrix(k, GATEID_I);
    }
    return prev;
}

/**
 * 
 * 
 * 
*/
MTBDD
mtbdd_stack_matrix(BDDVAR k, gate_id_t gateid) // What does low?
{
    // This function effectively does a Kronecker product gate \tensor below
    BDDVAR s, t;
    MTBDD u00, u01, u10, u11, low, high, res;

    // Even + uneven variable are used to encode the 4 values
    s = 2*k;
    t = s + 1;

    // Matrix U = [u00 u01
    //             u10 u11] encoded in a small tree
    u00 = mtbdd_makeleaf(MPC_TYPE, (size_t)gates[gateid][0]);
    u10 = mtbdd_makeleaf(MPC_TYPE, (size_t)gates[gateid][2]);
    u01 = mtbdd_makeleaf(MPC_TYPE, (size_t)gates[gateid][1]);
    u11 = mtbdd_makeleaf(MPC_TYPE, (size_t)gates[gateid][3]);

    low  = mtbdd_makenode(t, u00, u10);
    high = mtbdd_makenode(t, u01, u11);
    res  = mtbdd_makenode(s, low, high);

    // Propagate common factor on previous root amp to new root amp
    // AMP new_root_amp = wgt_mul(AADD_WEIGHT(below), AADD_WEIGHT(res));
    // res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}
