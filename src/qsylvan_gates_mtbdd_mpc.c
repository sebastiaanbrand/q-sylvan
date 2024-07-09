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
#include <sylvan_mpc.h>

/**
 * Definition of Quantum mpc values to fill the fixed and dynamic gates.
 * 
 * They are defined to be used globally to improve computation performance for the dynamic gates.
 * 
 * After gate initialization and use, call gate exit function to clear these mpc values.
 *  
*/

//
// For the fixed gates
//

// static MTBDD I_dd = MTBDD_ZERO; // moved to header file
// static MTBDD X_dd = MTBDD_ZERO;
// static MTBDD Y_dd = MTBDD_ZERO;
// static MTBDD Z_dd = MTBDD_ZERO;

// static MTBDD H_dd = MTBDD_ZERO;
// static MTBDD S_dd = MTBDD_ZERO;
// static MTBDD S_dag_dd = MTBDD_ZERO;
// static MTBDD T_dd = MTBDD_ZERO;
// static MTBDD T_dag_dd = MTBDD_ZERO;

// static MTBDD sqrt_X_dd = MTBDD_ZERO;
// static MTBDD sqrt_X_dag_dd = MTBDD_ZERO;
// static MTBDD sqrt_Y_dd = MTBDD_ZERO;
// static MTBDD sqrt_Y_dag_dd = MTBDD_ZERO;

static mpc_ptr mpc_pi;
static mpc_ptr mpc_zero;                      //  0.0 + i 0.0
static mpc_ptr mpc_re_one;                    //  1.0 + i 0.0
static mpc_ptr mpc_re_one_min;                //  0.0 - i 1.0
static mpc_ptr mpc_im_one;                    //  0.0 + i 1.0
static mpc_ptr mpc_im_one_min;                //  0.0 - i 0.0
static mpc_ptr mpc_half_half;                 //  0.5 + i 0.5
static mpc_ptr mpc_half_half_min;             //  0.5 - i 0.5
static mpc_ptr mpc_half_min_half;             // -0.5 + i 0.5
static mpc_ptr mpc_half_min_half_min;         // -0.5 - i 0.5

static mpc_ptr mpc_sqrt_2;                    //  sqrt(2)
static mpc_ptr mpc_res_sqrt_2;                //  1/sqrt(2)
static mpc_ptr mpc_res_sqrt_2_min;            // -1/sqrt(2)
static mpc_ptr mpc_res_sqrt_2_res_sqrt_2;     //  cos(pi/4) + i sin(pi/4) = 1/sqrt(2) + i/sqrt(2)
static mpc_ptr mpc_res_sqrt_2_res_sqrt_2_min; //  cos(pi/4) - i sin(pi/4) = 1/sqrt(2) - i/sqrt(2)

//
// For the dynamic initialized gates with theta Rx, Ry, Rz
//

#define MAX_QUBITS 128

static MTBDD R_dd[MAX_QUBITS], R_dag_dd[MAX_QUBITS];

static mpfr_t mpfr_pi;
static mpfr_t mpfr_zero;
static mpfr_t mpfr_theta;
static mpfr_t mpfr_theta_2;
static mpfr_t mpfr_min_theta_2;
static mpfr_t mpfr_cos_theta_2;
static mpfr_t mpfr_sin_theta_2;

static mpc_t mpc_cos_theta_2;           // cos(theta/2) + i 0.0
static mpc_t mpc_zero_sin_min_theta_2;  // 0.0 - i sin(theta/2)

static mpc_t mpc_sin_min_theta_2;       // -sin(theta/2) + i 0.0
static mpc_t mpc_sin_theta_2;           // sin(theta/2) + i 0.0

static mpfr_t mpfr_cos_min_theta_2;
static mpfr_t mpfr_sin_min_theta_2;
static mpc_t mpc_exp_min_theta_2;
static mpc_t mpc_exp_theta_2;

static mpfr_t mpfr_cos_theta;           // cos(theta)
static mpfr_t mpfr_sin_theta;           // sin(theta/2)
    
static mpc_t mpc_exp_theta;             // cos(theta) + i sin(theta)

static mpfr_t mpfr_phi;
static mpfr_t mpfr_lambda;
static mpfr_t mpfr_gam;

static mpfr_t mpfr_min_sin_theta_2, mpfr_sin_theta_2; // -sin(theta/2)

static mpc_t mpc_exp_lambda; 
static mpfr_t mpfr_cos_lambda, mpfr_sin_lambda; // cos(lambda) + i sin(lambda)
static mpc_t mpc_exp_lambda_mul_min_sin_theta_2, mpc_exp_lambda;
static mpfr_t mpfr_min_sin_theta_2;

static mpc_t mpc_exp_phi;
static mpfr_t mpfr_cos_phi, mpfr_sin_phi; // cos(phi) + i sin(phi)
static mpc_t mpc_exp_phi_mul_sin_theta_2, mpc_exp_phi;
static mpfr_t mpfr_sin_theta_2;

static mpc_t mpc_exp_gam; 
static mpfr_t mpfr_cos_gam, mpfr_sin_gam; // cos(gamma) + i sin(gamma)
static mpc_t mpc_exp_gam_mul_cos_theta_2, mpc_exp_gam;
static mpfr_t mpfr_cos_theta_2;


// All fixed 2x2 gates
//uint64_t G[nr_predef_gates+256+256][2][2]; // G[gateid][row][col]

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
    mpc_clear(mpc_pi);

    mpc_clear(mpc_pi);
    mpc_clear(mpc_zero);
    mpc_clear(mpc_re_one);
    mpc_clear(mpc_re_one_min);
    mpc_clear(mpc_im_one);
    mpc_clear(mpc_im_one_min);
    mpc_clear(mpc_half_half);
    mpc_clear(mpc_half_half_min);
    mpc_clear(mpc_half_min_half);
    mpc_clear(mpc_half_min_half_min);

    mpc_clear(mpc_sqrt_2);
    mpc_clear(mpc_res_sqrt_2);
    mpc_clear(mpc_res_sqrt_2_min);
    mpc_clear(mpc_res_sqrt_2_res_sqrt_2);
    mpc_clear(mpc_res_sqrt_2_res_sqrt_2_min);

    mpfr_clear(mpfr_theta);
    mpfr_clear(mpfr_theta_2);
    mpfr_clear(mpfr_cos_theta_2);
    mpfr_clear(mpfr_sin_theta_2);
    mpfr_clear(mpfr_cos_min_theta_2);
    mpfr_clear(mpfr_sin_min_theta_2);

    mpc_clear(mpc_exp_min_theta_2);
    mpc_clear(mpc_exp_theta_2);

    return;
}

/**
 * Initialize mpc variables to be used globally.
 */
void
mtbdd_gate_init_fixed_variables()
{
    mpfr_init2(mpfr_pi, MPC_PRECISION);
    mpfr_const_pi(mpfr_pi, MPC_ROUNDING);

    mpc_init2(mpc_pi, MPC_PRECISION);
    mpc_set_fr(mpc_pi, mpfr_pi, MPC_ROUNDING);

    mpc_assign(mpc_zero, 0.0, 0.0); // 0.0 + i 0.0
    mpc_assign(mpc_re_one, 1.0, 0.0); // 1.0 + i 0.0
    mpc_assign(mpc_re_one_min, -1.0, 0.0); // ...
    mpc_assign(mpc_im_one, 0.0, 1.0);
    mpc_assign(mpc_im_one_min, 0.0, -1.0);
    mpc_assign(mpc_half_half, 0.5, 0.5);
    mpc_assign(mpc_half_half_min, 0.5, -0.5);
    mpc_assign(mpc_half_min_half, -0.5, 0.5);
    mpc_assign(mpc_half_min_half_min, -0.5, -0.5);

    mpc_sqrt_assign(mpc_sqrt_2, 2.0, 0.0); // sqrt(2.0)

    mpc_init2(mpc_res_sqrt_2, MPC_PRECISION);
    mpc_div(mpc_res_sqrt_2, mpc_re_one, mpc_sqrt_2, MPC_ROUNDING); // 1/sqrt(2.0)

    mpc_init2(mpc_res_sqrt_2, MPC_PRECISION);
    mpc_div(mpc_res_sqrt_2_min, mpc_re_one_min, mpc_sqrt_2, MPC_ROUNDING); // -1/sqrt(2.0)

    mpc_init2(mpc_res_sqrt_2_res_sqrt_2, MPC_PRECISION);
    mpc_add(mpc_res_sqrt_2_res_sqrt_2, mpc_res_sqrt_2, mpc_res_sqrt_2, MPC_ROUNDING); // 1/sqrt(2.0) + i 1/sqrt(2.0)

    mpc_init2(mpc_res_sqrt_2_res_sqrt_2_min, MPC_PRECISION);
    mpc_add(mpc_res_sqrt_2_res_sqrt_2_min, mpc_res_sqrt_2, mpc_res_sqrt_2_min, MPC_ROUNDING); // 1/sqrt(2.0) - i 1/sqrt(2.0)

    return;
}

void
mtbdd_gate_init_dynamic_variables()
{
    mpfr_init2(mpfr_theta, MPC_PRECISION);
    mpfr_init2(mpfr_theta_2, MPC_PRECISION);
    mpfr_init2(mpfr_cos_theta_2, MPC_PRECISION);
    mpfr_init2(mpfr_sin_theta_2, MPC_PRECISION);
    mpfr_init2(mpfr_cos_min_theta_2, MPC_PRECISION);
    mpfr_init2(mpfr_sin_min_theta_2, MPC_PRECISION);

    mpc_init2(mpc_exp_min_theta_2, MPC_PRECISION);
    mpc_init2(mpc_exp_theta_2, MPC_PRECISION);

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

    G_arr[0][0] = mpc_re_one; 
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;   
    G_arr[1][1] = mpc_re_one;

    I_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_zero;   
    G_arr[0][1] = mpc_re_one;
    G_arr[1][0] = mpc_re_one; 
    G_arr[1][1] = mpc_zero;

    X_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_zero;   
    G_arr[0][1] = mpc_im_one_min;
    G_arr[1][0] = mpc_im_one; 
    G_arr[1][1] = mpc_zero;

    Y_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_re_one; 
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;   
    G_arr[1][1] = mpc_re_one_min;

    Z_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_res_sqrt_2; 
    G_arr[0][1] = mpc_res_sqrt_2;
    G_arr[1][0] = mpc_res_sqrt_2; 
    G_arr[1][1] = mpc_res_sqrt_2_min;

    H_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_re_one; 
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;   
    G_arr[1][1] = mpc_im_one;

    S_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_re_one; 
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;   
    G_arr[1][1] = mpc_im_one_min;

    S_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_re_one;  
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;    
    G_arr[1][1] = mpc_res_sqrt_2_res_sqrt_2;

    T_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_re_one;  
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;    
    G_arr[1][1] = mpc_res_sqrt_2_res_sqrt_2_min;

    T_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_half_half;     
    G_arr[0][1] = mpc_half_half_min;
    G_arr[1][0] = mpc_half_half_min; 
    G_arr[1][1] = mpc_half_half;

    sqrt_X_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_half_half_min; 
    G_arr[0][1] = mpc_half_half;
    G_arr[1][0] = mpc_half_half;     
    G_arr[1][1] = mpc_half_half_min;

    sqrt_X_dag_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_half_half;     
    G_arr[0][1] = mpc_half_min_half_min;
    G_arr[1][0] = mpc_half_half;     
    G_arr[1][1] = mpc_half_half;

    sqrt_Y_dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    G_arr[0][0] = mpc_half_half_min;  
    G_arr[0][1] = mpc_half_half_min; 
    G_arr[1][0] = mpc_half_min_half;  
    G_arr[1][1] = mpc_half_half_min; 

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
 *  Parameterized rotation gates Rx(theta), Ry(theta), Rz(theta)
 */

MTBDD
mtbdd_Rx(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(mpfr_zero, 0.0, MPC_ROUNDING);
    mpfr_set_d(mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(mpfr_theta_2, mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(mpfr_min_theta_2, mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(mpfr_cos_theta_2, mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(mpfr_sin_theta_2, mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(mpfr_cos_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(mpfr_sin_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(mpc_cos_theta_2, mpfr_cos_theta_2, mpfr_zero, MPC_ROUNDING); // cos(theta/2) + i 0.0
    mpc_set_fr_fr(mpc_zero_sin_min_theta_2, mpfr_zero, mpfr_sin_min_theta_2, MPC_ROUNDING); // 0.0 - i sin(theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = mpc_cos_theta_2;
    G_arr[0][1] = mpc_zero_sin_min_theta_2;
    G_arr[1][0] = mpc_zero_sin_min_theta_2;
    G_arr[1][1] = mpc_cos_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Ry(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(mpfr_theta_2, mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(mpfr_min_theta_2, mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(mpfr_cos_theta_2, mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(mpfr_sin_theta_2, mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(mpfr_cos_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(mpfr_sin_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(mpc_cos_theta_2, mpfr_cos_theta_2, mpfr_zero, MPC_ROUNDING); // cos(theta/2) + i 0.0
    mpc_set_fr_fr(mpc_sin_min_theta_2, mpfr_sin_min_theta_2, mpfr_zero, MPC_ROUNDING); // -sin(theta/2) + i 0.0
    mpc_set_fr_fr(mpc_sin_theta_2, mpfr_sin_theta_2, mpfr_zero, MPC_ROUNDING); // sin(theta/2) + i 0.0

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = mpc_cos_theta_2;
    G_arr[0][1] = mpc_sin_min_theta_2;
    G_arr[1][0] = mpc_sin_theta_2;
    G_arr[1][1] = mpc_cos_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Rz(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(mpfr_theta, theta, MPC_ROUNDING);
    mpfr_div_ui(mpfr_theta_2, mpfr_theta, 2, MPC_ROUNDING);  //  theta / 2
    mpfr_neg(mpfr_min_theta_2, mpfr_theta_2, MPC_ROUNDING);  // -theta / 2

    mpfr_cos(mpfr_cos_theta_2, mpfr_theta_2, MPC_ROUNDING);  // cos(theta/2)
    mpfr_sin(mpfr_sin_theta_2, mpfr_theta_2, MPC_ROUNDING);  // sin(theta/2)
    mpfr_cos(mpfr_cos_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // cos(-theta/2)
    mpfr_sin(mpfr_sin_min_theta_2, mpfr_min_theta_2, MPC_ROUNDING);  // sin(-theta/2)

    mpc_set_fr_fr(mpc_exp_theta_2, mpfr_cos_theta_2, mpfr_sin_theta_2, MPC_ROUNDING); // cos(theta/2) + i sin(theta/2)
    mpc_set_fr_fr(mpc_exp_min_theta_2, mpfr_cos_min_theta_2, mpfr_sin_min_theta_2, MPC_ROUNDING); // cos(-theta/2) + i sin(-theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = mpc_exp_min_theta_2;
    G_arr[0][1] = mpc_zero;
    G_arr[1][0] = mpc_zero;
    G_arr[1][1] = mpc_exp_theta_2;

    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_Phase(double theta) // TODO: change in mpfr!
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(mpfr_theta, theta, MPC_ROUNDING);

    mpfr_cos(mpfr_cos_theta, mpfr_theta, MPC_ROUNDING);  // cos(theta)
    mpfr_sin(mpfr_sin_theta, mpfr_theta, MPC_ROUNDING);  // sin(theta)

    mpc_set_fr_fr(mpc_exp_theta, mpfr_cos_theta, mpfr_sin_theta, MPC_ROUNDING); // cos(theta) + i sin(theta)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = mpc_re_one; // 1.0
    G_arr[0][1] = mpc_zero; // 0.0
    G_arr[1][0] = mpc_zero; // 0.0
    G_arr[1][1] = mpc_exp_theta; // exp(i theta) = cos(theta) + i sin(theta)
 
    dd = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

    free_matrix_array_mpc(G_arr, n);

    return dd;
}

MTBDD
mtbdd_U(double theta, double phi, double lambda)
{
    // Calculate complex numbers for gate elements

    mpfr_set_d(mpfr_theta, theta, MPC_ROUNDING);
    mpfr_set_d(mpfr_phi, phi, MPC_ROUNDING);
    mpfr_set_d(mpfr_lambda, lambda, MPC_ROUNDING);
    mpfr_set_d(mpfr_gam, phi + lambda, MPC_ROUNDING);

    mpfr_div_ui(mpfr_theta_2, mpfr_theta, 2, MPC_ROUNDING);           // theta / 2

    mpfr_cos(mpfr_cos_theta_2, mpfr_theta_2, MPC_ROUNDING);           // cos(theta/2)
    mpfr_sin(mpfr_sin_theta_2, mpfr_theta_2, MPC_ROUNDING);           // sin(theta/2)
    mpfr_neg(mpfr_min_sin_theta_2, mpfr_sin_theta_2, MPC_ROUNDING);   // -sin(theta/2)

    // phi cos / sin
    mpfr_cos(mpfr_cos_phi, mpfr_phi, MPC_ROUNDING);                   // cos(phi)
    mpfr_sin(mpfr_sin_phi, mpfr_phi, MPC_ROUNDING);                   // sin(phi)

    // lambda cos / sin
    mpfr_cos(mpfr_cos_lambda, mpfr_lambda, MPC_ROUNDING);             // cos(lambda)
    mpfr_sin(mpfr_sin_lambda, mpfr_lambda, MPC_ROUNDING);             // sin(lambda)

    // lambda cos / sin
    mpfr_cos(mpfr_cos_gam, mpfr_gam, MPC_ROUNDING);                   // cos(gamma)
    mpfr_sin(mpfr_sin_gam, mpfr_gam, MPC_ROUNDING);                   // sin(gamma)


    mpc_set_fr_fr(mpc_exp_phi, mpfr_cos_phi, mpfr_sin_phi, MPC_ROUNDING);                 // cos(phi) + i sin(phi)
    mpc_mul_fr(mpc_exp_phi_mul_sin_theta_2, mpc_exp_phi, mpfr_sin_theta_2, MPC_ROUNDING); // (cos(phi) + i sin(phi)) x sin(theta/2)

    mpc_set_fr_fr(mpc_exp_lambda, mpfr_cos_lambda, mpfr_sin_lambda, MPC_ROUNDING);                      // cos(lambda) + i sin(lambda)
    mpc_mul_fr(mpc_exp_lambda_mul_min_sin_theta_2, mpc_exp_lambda, mpfr_min_sin_theta_2, MPC_ROUNDING); // (cos(lambda) + i sin(lambda)) x -sin(theta/2)

    mpc_set_fr_fr(mpc_exp_gam, mpfr_cos_gam, mpfr_sin_gam, MPC_ROUNDING);                 // cos(gamma) + i sin(gamma)
    mpc_mul_fr(mpc_exp_gam_mul_cos_theta_2, mpc_exp_gam, mpfr_cos_theta_2, MPC_ROUNDING); // (cos(gamma) + i sin(gamma)) x cos(theta/2)

    // Convert to MTBDD

    MTBDD dd = MTBDD_ZERO;
    int n = 1;
    mpc_ptr **G_arr = NULL;
    allocate_matrix_array_mpc(&G_arr, n);

    G_arr[0][0] = mpc_cos_theta_2;                       // cos(theta/2)
    G_arr[0][1] = mpc_exp_lambda_mul_min_sin_theta_2;    // (cos(lambda) + i sin(lambda)) x -sin(theta/2)
    G_arr[1][0] = mpc_exp_phi_mul_sin_theta_2;           // (cos(phi) + i sin(phi)) x sin(theta/2)
    G_arr[1][1] = mpc_exp_gam_mul_cos_theta_2;           // (cos(phi+lambda) + i sin(phi+lambda)) x cos(theta/2)

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

    mpfr_mul_si(mpfr_2_pi, mpfr_pi, 2, MPC_ROUNDING); // 2 * pi
    mpfr_set_si(mpfr_2, 2, MPC_ROUNDING);

    // MTBDD R_dd[n], R_dag_dd[n];

    for (int k=0; k<=n; k++) {

        mpfr_pow_si(mpfr_2_power_k, mpfr_2, (long)k, MPC_ROUNDING); // 2 ^ k
        mpfr_div(mpfr_2_pi_div_2_power_k, mpfr_2_pi, mpfr_2_power_k, MPC_ROUNDING); // 2 * pi / (2 ^ k)
        mpfr_neg(mpfr_2_pi_div_2_power_k_min, mpfr_2_pi_div_2_power_k, MPC_ROUNDING); // -2 * pi / (2 ^ k)

        mpfr_cos(mpfr_cos_theta, mpfr_2_pi_div_2_power_k, MPC_ROUNDING);  // cos(theta)
        mpfr_sin(mpfr_sin_theta, mpfr_2_pi_div_2_power_k, MPC_ROUNDING);  // sin(theta)

        mpc_set_fr_fr(mpc_exp_theta, mpfr_cos_theta, mpfr_sin_theta, MPC_ROUNDING); // cos(theta) + i sin(theta)

        mpfr_cos(mpfr_cos_theta_min, mpfr_2_pi_div_2_power_k_min, MPC_ROUNDING);  // cos(-theta)
        mpfr_sin(mpfr_sin_theta_min, mpfr_2_pi_div_2_power_k_min, MPC_ROUNDING);  // sin(-theta)

        mpc_set_fr_fr(mpc_exp_theta_min, mpfr_cos_theta_min, mpfr_sin_theta_min, MPC_ROUNDING); // cos(-theta) + i sin(-theta)

        R_dd[k] = MTBDD_ZERO;
        R_dag_dd[k] = MTBDD_ZERO;

        mpc_ptr **G_arr = NULL;
        allocate_matrix_array_mpc(&G_arr, n);

        G_arr[0][0] = mpc_re_one;
        G_arr[0][1] = mpc_zero;
        G_arr[1][0] = mpc_zero;
        G_arr[1][1] = mpc_exp_theta;
 
        R_dd[k] = matrix_array_to_mtbdd_mpc(G_arr, n, ALTERNATE_ROW_FIRST_WISE_MODE);

        G_arr[0][0] = mpc_re_one;
        G_arr[0][1] = mpc_zero;
        G_arr[1][0] = mpc_zero;
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
