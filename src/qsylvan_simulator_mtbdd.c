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
mtbdd_create_all_zero_state(BDDVAR n)
{
    bool x[n];
    for (BDDVAR k=0; k<n; k++) x[k] = 0;
    return mtbdd_create_basis_state(n, x);
}

/**
 * Convert one state vector x[] = {0.0, ..., 1.0, 0.0} to MTBDD 
 * 
 * This results in an MTBDD with mpc complex leafs, deepness n
 * 
 *              x0
 *          x1      0.0+i0.0
 *  ...
 *          x(n-1)
 *  0.0+i0.0      1.0+i0.0
 * 
 */
MTBDD
mtbdd_create_basis_state(BDDVAR n, bool* x)
{
    MTBDD low, high, prev;

    mpc_t mpc_zero;
    mpc_init2(mpc_zero, MPC_PRECISION);
    mpc_assign(mpc_zero, 0.0, 0.0);

    mpc_t mpc_one;
    mpc_init2(mpc_one, MPC_PRECISION);
    mpc_assign(mpc_one, 1.0, 0.0);

    for (uint32_t k = n-1; k != (uint32_t)0; k--) {
        if (x[k] == 0) {
            if(k == (n-1))
                low = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_one);
            else
                low = prev;
            high = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_zero);
        }
        else if (x[k] == 1) {
            if(k == (n-1))
                low = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_zero);
            else
                low = prev;
            high = mtbdd_makeleaf(MPC_TYPE, (size_t)mpc_one);
        }
        prev = mtbdd_makenode(k, low, high);
    }

    mpc_clear(mpc_zero);
    mpc_clear(mpc_one);

    return prev;
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
