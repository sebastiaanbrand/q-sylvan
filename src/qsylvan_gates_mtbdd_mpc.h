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

#ifndef SYLVAN_GATES_MTBDD_MPC_H
#define SYLVAN_GATES_MTBDD_MPC_H

#include <stdint.h>
#include <sylvan_int.h>
#include <sylvan_mpc.h>

/**
 * 
 * Global defined mtbdd's of fixed gates.
 * 
 */

extern MTBDD I_dd;
extern MTBDD X_dd;
extern MTBDD Y_dd;
extern MTBDD Z_dd;

extern MTBDD H_dd;
extern MTBDD S_dd;
extern MTBDD S_dag_dd;
extern MTBDD T_dd;
extern MTBDD T_dag_dd;

extern MTBDD sqrt_X_dd;
extern MTBDD sqrt_X_dag_dd;
extern MTBDD sqrt_Y_dd;
extern MTBDD sqrt_Y_dag_dd;

extern MTBDD V00_dd;
extern MTBDD V11_dd;

struct mpc_variables_t {

    // Vars for fixed gates

    mpfr_t mpfr_zero;
    mpfr_t mpfr_pi;

    mpc_t mpc_zero;                      //  0.0 + i 0.0
    mpc_t mpc_re_one;                    //  1.0 + i 0.0
    mpc_t mpc_re_one_min;                // -1.0 + i 0.0
    mpc_t mpc_im_one;                    //  0.0 + i 1.0
    mpc_t mpc_im_one_min;                //  0.0 - i 1.0
    
    mpc_t mpc_half_half;                 //  0.5 + i 0.5
    mpc_t mpc_half_half_min;             //  0.5 - i 0.5
    mpc_t mpc_half_min_half;             // -0.5 + i 0.5
    mpc_t mpc_half_min_half_min;         // -0.5 - i 0.5

    mpc_t mpc_sqrt_2;                    //  sqrt(2)
    mpc_t mpc_res_sqrt_2;                //  1/sqrt(2)
    mpc_t mpc_im_res_sqrt_2;             //  i 1/sqrt(2)
    mpc_t mpc_res_sqrt_2_min;            // -1/sqrt(2)
    mpc_t mpc_im_res_sqrt_2_min;         // -i 1/sqrt(2)
    mpc_t mpc_res_sqrt_2_res_sqrt_2;     //  cos(pi/4) + i sin(pi/4) = 1/sqrt(2) + i/sqrt(2)
    mpc_t mpc_res_sqrt_2_res_sqrt_2_min; //  cos(pi/4) - i sin(pi/4) = 1/sqrt(2) - i/sqrt(2)


    // Vars for dynamic gates (parameterized gates)

    mpfr_t mpfr_theta;                   //  theta
    mpfr_t mpfr_theta_2;                 //  theta / 2
    mpfr_t mpfr_min_theta_2;             // -theta / 2
    mpfr_t mpfr_cos_theta_2;             //  cos(theta/2)
    mpfr_t mpfr_sin_theta_2;             //  sin(theta/2)

    mpc_t mpc_cos_theta_2;               //  cos(theta/2) + i 0.0
    mpc_t mpc_zero_sin_min_theta_2;      //  0.0 - i sin(theta/2)

    mpc_t mpc_sin_min_theta_2;           // -sin(theta/2) + i 0.0
    mpc_t mpc_sin_theta_2;               //  sin(theta/2) + i 0.0

    mpfr_t mpfr_cos_min_theta_2;         //  cos(-theta/2)
    mpfr_t mpfr_sin_min_theta_2;         //  sin(-theta/2)
    mpc_t mpc_exp_min_theta_2;           //  cos(-theta/2) + i sin(-theta/2)
    mpc_t mpc_exp_theta_2;               //  cos(theta/2) + i sin(theta/2)

    mpfr_t mpfr_cos_theta;               //  cos(theta)
    mpfr_t mpfr_sin_theta;               //  sin(theta)
    mpc_t mpc_exp_theta;                 //  cos(theta) + i sin(theta)

    mpfr_t mpfr_phi;                     //  phi
    mpfr_t mpfr_lambda;                  //  lambda
    mpfr_t mpfr_gam;                     //  gamma

    mpfr_t mpfr_min_sin_theta_2;         // -sin(theta/2)

    mpfr_t mpfr_cos_lambda;              //  cos(lambda)
    mpfr_t mpfr_sin_lambda;              //  sin(lambda)
    mpc_t mpc_exp_lambda;                //  cos(lambda) + i sin(lambda) 

    mpfr_t mpfr_cos_phi;                 //  cos(phi)
    mpfr_t mpfr_sin_phi;                 //  sin(phi)
    mpc_t mpc_exp_phi;                   //  cos(phi) + i sin(phi)

    mpfr_t mpfr_cos_gam;                 //  cos(gamma)
    mpfr_t mpfr_sin_gam;                 //  sin(gamma)
    mpc_t mpc_exp_gam;                   //  cos(gamma) + i sin(gamma) 

    mpc_t mpc_exp_phi_mul_sin_theta_2;          //  (cos(phi) + i sin(phi)) x cos(theta/2)
    mpc_t mpc_exp_lambda_mul_min_sin_theta_2;   //  (cos(lambda) + i sin(lambda)) x cos(theta/2)
    mpc_t mpc_exp_gam_mul_cos_theta_2;          //  (cos(gamma) + i sin(gamma)) x cos(theta/2)

};
extern struct mpc_variables_t g;

/**
 * 
 * Functions
 * 
 */

// "Constructor"
void 
mtbdd_gates_init_mpc();

// "Destructor"
void
mtbdd_gate_exit_mpc();

//
// TODO: optimize computation speed by re-use of dynamic gate
//
// The reason why these are initialized beforehand instead of on-demand is that 
// we would like a (for example) pi/16 gate to always have the same unique ID 
// throughout the entire run of the circuit.
//

void
mtbdd_gate_init_fixed_variables();

void
mtbdd_gate_init_dynamic_variables();

void
mtbdd_fixed_gates_init_mpc();

void 
mtbdd_phase_gates_init(int n);

/**
 * Rotation around x-axis with angle theta.
 */
MTBDD 
mtbdd_Rx(double theta);

/**
 * Rotation around y-axis with angle theta.
 */
MTBDD 
mtbdd_Ry(double theta);

/**
 * Rotation around z-axis with angle theta.
 */
MTBDD 
mtbdd_Rz(double theta);

/**
 * Rotation around z-axis with angle theta (but different global phase than Rz)
 */
MTBDD 
mtbdd_Phase(double theta);

/**
 * Generic single-qubit rotation gate with 3 Euler angles.
 */
MTBDD 
mtbdd_U(double theta, double phi, double lambda);

#endif
