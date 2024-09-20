/**
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

#include <qsylvan_gates_mtbdd_mpc.h>
//#include <sylvan_int.h>
//#include <sylvan_edge_weights_complex.h>


/**
 * Definition of Quantum mpc values to fill the fixed and dynamic gates.
 * 
 * They are defined to be used globally to improve computation performance for the dynamic gates.
 * 
 * After gate initialization and use, call gate exit function to clear these mpc values.
 *  
 */

#define MAX_QUBITS 128

static MTBDD R_dd[MAX_QUBITS], R_dag_dd[MAX_QUBITS];

// All fixed 2x2 gates
// uint64_t G[nr_predef_gates+256+256][2][2]; // G[gateid][row][col]

/**
 * 
 * Global defined mtbdd's of fixed gates.
 * 
 */
MTBDD I_dd = MTBDD_ZERO;
MTBDD X_dd = MTBDD_ZERO;
MTBDD Y_dd = MTBDD_ZERO;
MTBDD Z_dd = MTBDD_ZERO;

MTBDD H_dd = MTBDD_ZERO;
MTBDD S_dd = MTBDD_ZERO;
MTBDD S_dag_dd = MTBDD_ZERO;
MTBDD T_dd = MTBDD_ZERO;
MTBDD T_dag_dd = MTBDD_ZERO;

MTBDD sqrt_X_dd = MTBDD_ZERO;
MTBDD sqrt_X_dag_dd = MTBDD_ZERO;
MTBDD sqrt_Y_dd = MTBDD_ZERO;
MTBDD sqrt_Y_dag_dd = MTBDD_ZERO;

/**
 * Collection of mpc variables to be used globally.
 */
struct mpc_variables_t g; // g = global

/**
 *  "Constructor"
 */
void
mtbdd_gates_init_mpc()
{
    mtbdd_gate_init_fixed_variables();
    mtbdd_gate_init_dynamic_variables();
    mtbdd_fixed_gates_init_mpc();

    return;
}

/**
 *  "Destructor"
 */
void
mtbdd_gate_exit_mpc()
{
    mpfr_clear(g.mpfr_pi);

    mpc_clear(g.mpc_zero);
    mpc_clear(g.mpc_re_one);
    mpc_clear(g.mpc_re_one_min);
    mpc_clear(g.mpc_im_one);
    mpc_clear(g.mpc_im_one_min);
    mpc_clear(g.mpc_half_half);
    mpc_clear(g.mpc_half_half_min);
    mpc_clear(g.mpc_half_min_half);
    mpc_clear(g.mpc_half_min_half_min);

    mpc_clear(g.mpc_sqrt_2);
    mpc_clear(g.mpc_res_sqrt_2);
    mpc_clear(g.mpc_res_sqrt_2_min);
    mpc_clear(g.mpc_res_sqrt_2_res_sqrt_2);
    mpc_clear(g.mpc_res_sqrt_2_res_sqrt_2_min);

    mpfr_clear(g.mpfr_theta);
    mpfr_clear(g.mpfr_theta_2);
    mpfr_clear(g.mpfr_cos_theta_2);
    mpfr_clear(g.mpfr_sin_theta_2);
    mpfr_clear(g.mpfr_cos_min_theta_2);
    mpfr_clear(g.mpfr_sin_min_theta_2);

    mpc_clear(g.mpc_exp_min_theta_2);
    mpc_clear(g.mpc_exp_theta_2);

    return;
}

/**
 * Declare several constants globally (g.) as mpc type 
 */
void
mtbdd_gate_init_fixed_variables()
{
    // pi, 0, 1, -1, ...
    mpfr_init2(g.mpfr_pi, MPC_PRECISION);
    mpfr_const_pi(g.mpfr_pi, MPC_ROUNDING);

    mpc_assign(g.mpc_zero, 0.0, 0.0);
    mpc_assign(g.mpc_re_one, 1.0, 0.0);
    mpc_assign(g.mpc_re_one_min, -1.0, 0.0);
    mpc_assign(g.mpc_im_one, 0.0, 1.0);
    mpc_assign(g.mpc_im_one_min, 0.0, -1.0);
    mpc_assign(g.mpc_half_half, 0.5, 0.5);
    mpc_assign(g.mpc_half_half_min, 0.5, -0.5);
    mpc_assign(g.mpc_half_min_half, -0.5, 0.5);
    mpc_assign(g.mpc_half_min_half_min, -0.5, -0.5);

    // sqrt(2.0)
    mpc_sqrt_assign(g.mpc_sqrt_2, 2.0, 0.0);

    // 1/sqrt(2.0)
    mpc_init2(g.mpc_res_sqrt_2, MPC_PRECISION);
    mpc_div(g.mpc_res_sqrt_2, g.mpc_re_one, g.mpc_sqrt_2, MPC_ROUNDING);

    // -1/sqrt(2.0)
    mpc_init2(g.mpc_res_sqrt_2_min, MPC_PRECISION);
    mpc_div(g.mpc_res_sqrt_2_min, g.mpc_re_one_min, g.mpc_sqrt_2, MPC_ROUNDING); 

    // 1/sqrt(2.0) + i 1/sqrt(2.0)
    mpc_init2(g.mpc_res_sqrt_2_res_sqrt_2, MPC_PRECISION);
    mpc_add(g.mpc_res_sqrt_2_res_sqrt_2, g.mpc_res_sqrt_2, g.mpc_res_sqrt_2, MPC_ROUNDING); 

    // 1/sqrt(2.0) - i 1/sqrt(2.0)
    mpc_init2(g.mpc_res_sqrt_2_res_sqrt_2_min, MPC_PRECISION);
    mpc_add(g.mpc_res_sqrt_2_res_sqrt_2_min, g.mpc_res_sqrt_2, g.mpc_res_sqrt_2_min, MPC_ROUNDING); 

    return;
}

/**
 * Initialize vars for dynamic composition of gates to be used globally (g.) as mpc type
 */
void
mtbdd_gate_init_dynamic_variables()
{
    mpfr_init2(g.mpfr_theta, MPC_PRECISION);
    mpfr_init2(g.mpfr_theta_2, MPC_PRECISION);
    mpfr_init2(g.mpfr_cos_theta_2, MPC_PRECISION);
    mpfr_init2(g.mpfr_sin_theta_2, MPC_PRECISION);
    mpfr_init2(g.mpfr_cos_min_theta_2, MPC_PRECISION);
    mpfr_init2(g.mpfr_sin_min_theta_2, MPC_PRECISION);

    mpc_init2(g.mpc_exp_min_theta_2, MPC_PRECISION);
    mpc_init2(g.mpc_exp_theta_2, MPC_PRECISION);

    return;
}

/**
 * Initialize all fixed gates with mpc type.
 */
void
mtbdd_fixed_gates_init_mpc()
{
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_re_one; 
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;   
    G_arr[1][1] = g.mpc_re_one;

    I_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_zero;   
    G_arr[0][1] = g.mpc_re_one;
    G_arr[1][0] = g.mpc_re_one; 
    G_arr[1][1] = g.mpc_zero;

    X_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_zero;   
    G_arr[0][1] = g.mpc_im_one_min;
    G_arr[1][0] = g.mpc_im_one; 
    G_arr[1][1] = g.mpc_zero;

    Y_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_re_one; 
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;   
    G_arr[1][1] = g.mpc_re_one_min;

    Z_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_res_sqrt_2; 
    G_arr[0][1] = g.mpc_res_sqrt_2;
    G_arr[1][0] = g.mpc_res_sqrt_2; 
    G_arr[1][1] = g.mpc_res_sqrt_2_min;

    H_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_re_one; 
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;   
    G_arr[1][1] = g.mpc_im_one;

    S_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_re_one; 
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;   
    G_arr[1][1] = g.mpc_im_one_min;

    S_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_re_one;  
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;    
    G_arr[1][1] = g.mpc_res_sqrt_2_res_sqrt_2;

    T_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_re_one;  
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;    
    G_arr[1][1] = g.mpc_res_sqrt_2_res_sqrt_2_min;

    T_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_half_half;     
    G_arr[0][1] = g.mpc_half_half_min;
    G_arr[1][0] = g.mpc_half_half_min; 
    G_arr[1][1] = g.mpc_half_half;

    sqrt_X_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_half_half_min; 
    G_arr[0][1] = g.mpc_half_half;
    G_arr[1][0] = g.mpc_half_half;     
    G_arr[1][1] = g.mpc_half_half_min;

    sqrt_X_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_half_half;     
    G_arr[0][1] = g.mpc_half_min_half_min;
    G_arr[1][0] = g.mpc_half_half;     
    G_arr[1][1] = g.mpc_half_half;

    sqrt_Y_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = g.mpc_half_half_min;  
    G_arr[0][1] = g.mpc_half_half_min; 
    G_arr[1][0] = g.mpc_half_min_half;  
    G_arr[1][1] = g.mpc_half_half_min; 

    sqrt_Y_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return;
}

/*
    k = GATEID_I;
    G[k][0][0] = mpc_re_one; G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;   G[k][1][1] = mpc_re_one;

    k = GATEID_X;
    G[k][0][0] = mpc_zero;   G[k][0][1] = mpc_re_one;
    G[k][1][0] = mpc_re_one; G[k][1][1] = mpc_zero;

    k = GATEID_Y;
    G[k][0][0] = mpc_zero;   G[k][0][1] = mpc_im_one_min;
    G[k][1][0] = mpc_im_one; G[k][1][1] = mpc_zero;

    k = GATEID_Z;
    G[k][0][0] = mpc_re_one; G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;   G[k][1][1] = mpc_re_one_min;

    k = GATEID_H;
    G[k][0][0] = mpc_res_sqrt_2; G[k][0][1] = mpc_res_sqrt_2;
    G[k][1][0] = mpc_res_sqrt_2; G[k][1][1] = mpc_res_sqrt_2_min;

    k = GATEID_S;
    G[k][0][0] = mpc_re_one; G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;   G[k][1][1] = mpc_im_one;

    k = GATEID_Sdag;
    G[k][0][0] = mpc_re_one; G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;   G[k][1][1] = mpc_im_one_min;

    k = GATEID_T;
    G[k][0][0] = mpc_re_one;  G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;    G[k][1][1] = mpc_res_sqrt_2_res_sqrt_2;

    k = GATEID_Tdag;
    G[k][0][0] = mpc_re_one;  G[k][0][1] = mpc_zero;
    G[k][1][0] = mpc_zero;    G[k][1][1] = mpc_res_sqrt_2_res_sqrt_2_min;

    k = GATEID_sqrtX;
    G[k][0][0] = mpc_half_half;     G[k][0][1] = mpc_half_half_min;
    G[k][1][0] = mpc_half_half_min; G[k][1][1] = mpc_half_half;

    k = GATEID_sqrtXdag;
    G[k][0][0] = mpc_half_half_min; G[k][0][1] = mpc_half_half;
    G[k][1][0] = mpc_half_half;     G[k][1][1] = mpc_half_half_min;

    k = GATEID_sqrtY;
    G[k][0][0] = mpc_half_half;     G[k][0][1] = mpc_half_min_half_min;
    G[k][1][0] = mpc_half_half;     G[k][1][1] = mpc_half_half;

    k = GATEID_sqrtYdag;
    G[k][0][0] = mpc_half_half_min;  G[k][0][1] = mpc_half_half_min; 
    G[k][1][0] = mpc_half_min_half;  G[k][1][1] = mpc_half_half_min; 
*/

/*
    qmdd_phase_gates_init(255);

    // init dynamic gate 
    // (necessary when qmdd_gates_init() is called after gc to re-init all gates)
    k = GATEID_dynamic;
    gates[k][0] = weight_lookup(&dynamic_gate[0]);
    gates[k][1] = weight_lookup(&dynamic_gate[1]);
    gates[k][2] = weight_lookup(&dynamic_gate[2]);
    gates[k][3] = weight_lookup(&dynamic_gate[3]);

    return;
}
*/

/**
 *  Parameterized dynamically computed rotation gates Rx(theta), Ry(theta), Rz(theta)
 */

MTBDD
mtbdd_Rx(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(g.mpfr_zero, 0.0, MPC_ROUNDING);
    mpfr_set_d(g.mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(g.mpfr_theta_2, g.mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(g.mpfr_min_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(g.mpfr_cos_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(g.mpfr_sin_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(g.mpfr_cos_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(g.mpfr_sin_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(g.mpc_cos_theta_2, g.mpfr_cos_theta_2, g.mpfr_zero, MPC_ROUNDING); // cos(theta/2) + i 0.0
    mpc_set_fr_fr(g.mpc_zero_sin_min_theta_2, g.mpfr_zero, g.mpfr_sin_min_theta_2, MPC_ROUNDING); // 0.0 - i sin(theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_cos_theta_2;
    G_arr[0][1] = g.mpc_zero_sin_min_theta_2;
    G_arr[1][0] = g.mpc_zero_sin_min_theta_2;
    G_arr[1][1] = g.mpc_cos_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Ry(double theta) // TODO: change to mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(g.mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(g.mpfr_theta_2, g.mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(g.mpfr_min_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(g.mpfr_cos_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(g.mpfr_sin_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(g.mpfr_cos_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(g.mpfr_sin_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(g.mpc_cos_theta_2, g.mpfr_cos_theta_2, g.mpfr_zero, MPC_ROUNDING); // cos(theta/2) + i 0.0
    mpc_set_fr_fr(g.mpc_sin_min_theta_2, g.mpfr_sin_min_theta_2, g.mpfr_zero, MPC_ROUNDING); // -sin(theta/2) + i 0.0
    mpc_set_fr_fr(g.mpc_sin_theta_2, g.mpfr_sin_theta_2, g.mpfr_zero, MPC_ROUNDING); // sin(theta/2) + i 0.0

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_cos_theta_2;
    G_arr[0][1] = g.mpc_sin_min_theta_2;
    G_arr[1][0] = g.mpc_sin_theta_2;
    G_arr[1][1] = g.mpc_cos_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Rz(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(g.mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(g.mpfr_theta_2, g.mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(g.mpfr_min_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(g.mpfr_cos_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(g.mpfr_sin_theta_2, g.mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(g.mpfr_cos_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(g.mpfr_sin_min_theta_2, g.mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(g.mpc_exp_theta_2, g.mpfr_cos_theta_2, g.mpfr_sin_theta_2, MPC_ROUNDING); // cos(theta/2) + i sin(theta/2)
    mpc_set_fr_fr(g.mpc_exp_min_theta_2, g.mpfr_cos_min_theta_2, g.mpfr_sin_min_theta_2, MPC_ROUNDING); // cos(-theta/2) + i sin(-theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_exp_min_theta_2;
    G_arr[0][1] = g.mpc_zero;
    G_arr[1][0] = g.mpc_zero;
    G_arr[1][1] = g.mpc_exp_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Phase(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(g.mpfr_theta, theta, MPC_ROUNDING);

    mpfr_cos(g.mpfr_cos_theta, g.mpfr_theta, MPC_ROUNDING);  // cos(theta)
    mpfr_sin(g.mpfr_sin_theta, g.mpfr_theta, MPC_ROUNDING);  // sin(theta)

    mpc_set_fr_fr(g.mpc_exp_theta, g.mpfr_cos_theta, g.mpfr_sin_theta, MPC_ROUNDING); // cos(theta) + i sin(theta)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_re_one; // 1.0
    G_arr[0][1] = g.mpc_zero; // 0.0
    G_arr[1][0] = g.mpc_zero; // 0.0
    G_arr[1][1] = g.mpc_exp_theta; // exp(i theta) = cos(theta) + i sin(theta)
 
    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_U(double theta, double phi, double lambda)
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(g.mpfr_theta, theta, MPC_ROUNDING);
    mpfr_set_d(g.mpfr_phi, phi, MPC_ROUNDING);
    mpfr_set_d(g.mpfr_lambda, lambda, MPC_ROUNDING);
    mpfr_set_d(g.mpfr_gam, phi + lambda, MPC_ROUNDING);

    mpfr_div_ui(g.mpfr_theta_2, g.mpfr_theta, 2, MPC_ROUNDING);           // theta / 2

    mpfr_cos(g.mpfr_cos_theta_2, g.mpfr_theta_2, MPC_ROUNDING);           // cos(theta/2)
    mpfr_sin(g.mpfr_sin_theta_2, g.mpfr_theta_2, MPC_ROUNDING);           // sin(theta/2)
    mpfr_neg(g.mpfr_min_sin_theta_2, g.mpfr_sin_theta_2, MPC_ROUNDING);   // -sin(theta/2)

    // phi cos / sin
    mpfr_cos(g.mpfr_cos_phi, g.mpfr_phi, MPC_ROUNDING);                   // cos(phi)
    mpfr_sin(g.mpfr_sin_phi, g.mpfr_phi, MPC_ROUNDING);                   // sin(phi)

    // lambda cos / sin
    mpfr_cos(g.mpfr_cos_lambda, g.mpfr_lambda, MPC_ROUNDING);             // cos(lambda)
    mpfr_sin(g.mpfr_sin_lambda, g.mpfr_lambda, MPC_ROUNDING);             // sin(lambda)

    // lambda cos / sin
    mpfr_cos(g.mpfr_cos_gam, g.mpfr_gam, MPC_ROUNDING);                   // cos(gamma)
    mpfr_sin(g.mpfr_sin_gam, g.mpfr_gam, MPC_ROUNDING);                   // sin(gamma)


    mpc_set_fr_fr(g.mpc_exp_phi, g.mpfr_cos_phi, g.mpfr_sin_phi, MPC_ROUNDING);                 // cos(phi) + i sin(phi)
    mpc_mul_fr(g.mpc_exp_phi_mul_sin_theta_2, g.mpc_exp_phi, g.mpfr_sin_theta_2, MPC_ROUNDING); // (cos(phi) + i sin(phi)) x sin(theta/2)

    mpc_set_fr_fr(g.mpc_exp_lambda, g.mpfr_cos_lambda, g.mpfr_sin_lambda, MPC_ROUNDING);                      // cos(lambda) + i sin(lambda)
    mpc_mul_fr(g.mpc_exp_lambda_mul_min_sin_theta_2, g.mpc_exp_lambda, g.mpfr_min_sin_theta_2, MPC_ROUNDING); // (cos(lambda) + i sin(lambda)) x -sin(theta/2)

    mpc_set_fr_fr(g.mpc_exp_gam, g.mpfr_cos_gam, g.mpfr_sin_gam, MPC_ROUNDING);                 // cos(gamma) + i sin(gamma)
    mpc_mul_fr(g.mpc_exp_gam_mul_cos_theta_2, g.mpc_exp_gam, g.mpfr_cos_theta_2, MPC_ROUNDING); // (cos(gamma) + i sin(gamma)) x cos(theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = g.mpc_cos_theta_2;                       // cos(theta/2)
    G_arr[0][1] = g.mpc_exp_lambda_mul_min_sin_theta_2;    // (cos(lambda) + i sin(lambda)) x -sin(theta/2)
    G_arr[1][0] = g.mpc_exp_phi_mul_sin_theta_2;           // (cos(phi) + i sin(phi)) x sin(theta/2)
    G_arr[1][1] = g.mpc_exp_gam_mul_cos_theta_2;           // (cos(phi+lambda) + i sin(phi+lambda)) x cos(theta/2)

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

/*
void
qmdd_phase_gates_init(int n) // n is number of qubits?
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    fl_t angle;
    complex_t cartesian;
*/

void
mtbdd_phase_gates_init(int n)
{
    mpfr_t mpfr_2_pi;
    mpfr_t mpfr_2;

    mpfr_t mpfr_2_power_k;
    mpfr_t mpfr_2_pi_div_2_power_k;
    mpfr_t mpfr_2_pi_div_2_power_k_min;

    mpfr_t mpfr_cos_theta_min;
    mpfr_t mpfr_sin_theta_min;

    mpc_t mpc_exp_theta_min;

    mpfr_init2(mpfr_2_pi, MPC_PRECISION);
    mpfr_init2(mpfr_2, MPC_PRECISION);
    mpfr_init2(mpfr_2_power_k, MPC_PRECISION);
    mpfr_init2(mpfr_2_pi_div_2_power_k, MPC_PRECISION);

    mpfr_mul_si(mpfr_2_pi, g.mpfr_pi, 2, MPC_ROUNDING); // 2 * pi
    mpfr_set_si(mpfr_2, 2, MPC_ROUNDING);

    // MTBDD R_dd[n], R_dag_dd[n];

    for (int k=0; k<=n; k++) {

        mpfr_pow_si(mpfr_2_power_k, mpfr_2, (long)k, MPC_ROUNDING); // 2 ^ k
        mpfr_div(mpfr_2_pi_div_2_power_k, mpfr_2_pi, mpfr_2_power_k, MPC_ROUNDING); // 2 * pi / (2 ^ k)
        mpfr_neg(mpfr_2_pi_div_2_power_k_min, mpfr_2_pi_div_2_power_k, MPC_ROUNDING); // -2 * pi / (2 ^ k)

        mpfr_cos(g.mpfr_cos_theta, mpfr_2_pi_div_2_power_k, MPC_ROUNDING);  // cos(theta)
        mpfr_sin(g.mpfr_sin_theta, mpfr_2_pi_div_2_power_k, MPC_ROUNDING);  // sin(theta)

        mpc_set_fr_fr(g.mpc_exp_theta, g.mpfr_cos_theta, g.mpfr_sin_theta, MPC_ROUNDING); // cos(theta) + i sin(theta)

        mpfr_cos(mpfr_cos_theta_min, mpfr_2_pi_div_2_power_k_min, MPC_ROUNDING);  // cos(-theta)
        mpfr_sin(mpfr_sin_theta_min, mpfr_2_pi_div_2_power_k_min, MPC_ROUNDING);  // sin(-theta)

        mpc_set_fr_fr(mpc_exp_theta_min, mpfr_cos_theta_min, mpfr_sin_theta_min, MPC_ROUNDING); // cos(-theta) + i sin(-theta)

        R_dd[k] = MTBDD_ZERO;
        R_dag_dd[k] = MTBDD_ZERO;

        mpc_ptr **G_arr = NULL;
        allocate_matrix_array_mpc(&G_arr, n);

        G_arr[0][0] = g.mpc_re_one;
        G_arr[0][1] = g.mpc_zero;
        G_arr[1][0] = g.mpc_zero;
        G_arr[1][1] = g.mpc_exp_theta;
 
        R_dd[k] = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

        G_arr[0][0] = g.mpc_re_one;
        G_arr[0][1] = g.mpc_zero;
        G_arr[1][0] = g.mpc_zero;
        G_arr[1][1] = mpc_exp_theta_min;
 
        R_dag_dd[k] = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

        free_matrix_array_mpc(G_arr, n);

    }

    mpfr_clear(mpfr_2_pi);
    mpfr_clear(mpfr_2);

    mpfr_clear(mpfr_2_power_k);
    mpfr_clear(mpfr_2_pi_div_2_power_k);
    mpfr_clear(mpfr_2_pi_div_2_power_k_min);

    mpfr_clear(mpfr_cos_theta_min);
    mpfr_clear(mpfr_sin_theta_min);

    mpc_clear(mpc_exp_theta_min);

/*
    for (int k=0; k<=n; k++) {

        // forward rotation
        angle = 2 * Pi / (1<<k); // (2 * pi) / (2 ^ k), so if n=2 (2 qubits), k = 0,1,2, angle = 2 * pi / (1,2,4) = 2pi, pi, pi/2

        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = AADD_ONE;  gates[gate_id][1] = AADD_ZERO;
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = AADD_ONE;  gates[gate_id][1] = AADD_ZERO;
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);
    }
*/

}
