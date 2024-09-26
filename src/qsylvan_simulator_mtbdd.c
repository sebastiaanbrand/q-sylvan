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
#include <qsylvan_gates_mtbdd_mpc.h>

/**
 * Convert zero state vector x[] = {0.0, 0.0, ..., 0.0} to MTBDD 
 * 
 * This results in an MTBDD with mpc complex zero leafs, deepness n
 * 
 *              x0
 *          x1      0.0+i0.0
 *  ...
 *          x(n-1)
 *
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
mtbdd_create_basis_state_double(BDDVAR n, bool* x) // TODO: convert double to complex_t
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
        //printf("x[%d] = %d, var = %d\n", i, x[i], var);

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
    //mpc_init2(mpc_zero, MPC_PRECISION);
    mpc_assign(mpc_zero, 0.0, 0.0);

    mpc_t mpc_one;
    //mpc_init2(mpc_one, MPC_PRECISION);
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
        //printf("x[%d] = %d, var = %d\n", i, x[i], var);

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
 * 
 * Create a single Qubit gate surrounded by I gates:
 * 
 *  q(n-1) ----- I(n-1) -----
 *  q(n-2) ----- I(n-2) -----
 *  q( . ) ----- I( . ) -----
 *  q( t ) ----- G( t ) -----
 *  q( . ) ----- I( . ) -----
 *  q( 0 ) ----- I( 0 ) -----
 * 
 *  q(i) = {0,1}, a qubit
 * 
 *  n is number of qubits
 *  t is index of qubit connected to the gate G
 * 
 *  (I(0) x ... x G(t) x ... x I(n-1) ) | q(n-1)q(n-2)...q(t)...q0>
 *
 *  a x b is the tensor of a and b, with a and b matrices / vectors.
 *  
 */
MTBDD
mtbdd_create_single_gate_for_qubits_mpc(BDDVAR n, BDDVAR t, MTBDD I_dd, MTBDD G_dd) //gate_id_t gateid)
{
    //MTBDD dd = I_dd;
    MTBDD dd = mtbdd_makeleaf(MPC_TYPE, (uint64_t)g.mpc_re_one);

    for(uint32_t k=0; k < n; k++)
    {
        if(k==t)
            dd = mtbdd_tensor_prod(G_dd, dd, 2); // Two vars added, so third argument = 2
        else
            dd = mtbdd_tensor_prod(I_dd, dd, 2); // Two vars added
    }

    return dd;
}

double
mtbdd_getnorm_mpc(MTBDD dd, size_t nvars) // L2 norm, in accordance with the satcount function in sylvan_mtbdd.c
{
    /* Trivial cases */
    //if (dd == mtbdd_false) 
    //    return 0.0;

    if (mtbdd_isleaf(dd)) {

        mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(dd), MPC_ROUNDING);

        mpc_ptr mpc_value = (mpc_ptr)mtbdd_getvalue(dd);
        
        //printf("%p\n", mpc_value);
        //mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, (mpc_ptr)mtbdd_getvalue(dd), MPC_ROUNDING);

        mpfr_t mpfr_abs_value;
        mpfr_init2(mpfr_abs_value, MPC_PRECISION);
        mpc_abs(mpfr_abs_value, mpc_value, MPC_ROUNDING);
    
        double double_value = mpfr_get_d(mpfr_abs_value, MPC_ROUNDING);
        mpfr_clear(mpfr_abs_value);

        //printf("getnorm_mpc = %f\n", double_value * double_value * pow(2.0, nvars));
        
        return double_value * double_value * pow(2.0, nvars);
        // return powl(2.0L, nvars);
    }

    /* Perhaps execute garbage collection */
    //sylvan_gc_test(); reads only the cache, can be removed

    union {   // copy bitvalues of double into 64 bit integer for cache
        double d;
        uint64_t s;
    } hack;

    /* Consult cache */
    if (cache_get3(CACHE_MTBDD_GETNORM_MPC, dd, 0, nvars, &hack.s)) {
        return hack.d;
    }

    double high = mtbdd_getnorm_mpc(mtbdd_gethigh(dd), nvars-1);
    double low = mtbdd_getnorm_mpc(mtbdd_getlow(dd), nvars-1);
    hack.d = low + high;

    cache_put3(CACHE_MTBDD_GETNORM_MPC, dd, 0, nvars, hack.s);

    return hack.d;
}




/*    
    int leaf_var_dd_I;

    MTBDD dd_I = mtbdd_of_identity_matrix(n);  // n is number of qubits, how do we know?
    MTBDD dd_G = mtbdd_of_gate_matrix(GATEID_X);
    MTBDD dd_gate = mtbdd_tensor_prod(dd_I, dd_G, leaf_var_dd_I);

    //MTBDD dd_state = mtbdd_create_state(state);

    return mtbdd_matvec_mult(dd_gate, state, (1 << n), 0);

    mpc_t mpc_zero;
    mpc_init2(mpc_zero, MPC_PRECISION);
    mpc_assign(mpc_zero, 0.0, 0.0);

    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    MTBDD node = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_one);

    for (int k = n-1; k >= 0; k--) {

        if ((unsigned int)k == t)
            node = mtbdd_stack_matrix(k, gateid);
        else
            node = mtbdd_stack_matrix(k, GATEID_I);
    }
    return node;
}

MTBDD
mtbdd_stack_matrix(BDDVAR k, gate_id_t gateid)
{
    // This function effectively does a Kronecker product gate \tensor below
    BDDVAR s, t;
    MTBDD u00, u01, u10, u11, low, high, res;

    // Even + uneven variable are used to encode the 4 values
    s = 2 * k;
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
*/

