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

/**
 * 
 * Global defined mtbdd's of fixed gates.
 * 
 */

static MTBDD I_dd = MTBDD_ZERO;
static MTBDD X_dd = MTBDD_ZERO;
static MTBDD Y_dd = MTBDD_ZERO;
static MTBDD Z_dd = MTBDD_ZERO;

static MTBDD H_dd = MTBDD_ZERO;
static MTBDD S_dd = MTBDD_ZERO;
static MTBDD S_dag_dd = MTBDD_ZERO;
static MTBDD T_dd = MTBDD_ZERO;
static MTBDD T_dag_dd = MTBDD_ZERO;

static MTBDD sqrt_X_dd = MTBDD_ZERO;
static MTBDD sqrt_X_dag_dd = MTBDD_ZERO;
static MTBDD sqrt_Y_dd = MTBDD_ZERO;
static MTBDD sqrt_Y_dag_dd = MTBDD_ZERO;

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

// The reason why these are initialized beforehand instead of on-demand is that 
// we would like a (for example) pi/16 gate to always have the same unique ID 
// throughout the entire run of the circuit.

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
