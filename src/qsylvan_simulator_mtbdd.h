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

#ifndef QSYLVAN_SIMULATOR_MTBDD_H
#define QSYLVAN_SIMULATOR_MTBDD_H

//#include <sylvan_int.h>
#include <qsylvan.h>
//#include <qsylvan_gates_mtbdd_mpc.h>
//#include <qsylvan_gates.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

//static const uint64_t MTBDD_TERMINAL = 1;
//static const BDDVAR   MTBDD_INVALID_VAR = UINT8_MAX;

/**
 * Creates an MTBDD for an n-qubit state |00...0>.
 * 
 * @param n Number of qubits.
 * 
 * @return An MTBDD encoding the n-qubit state |00..0>.
 */
MTBDD mtbdd_create_all_zero_state_double(BDDVAR n); //TODO: extend with complex_t
MTBDD mtbdd_create_all_zero_state_mpc(BDDVAR n);

/**
 * Creates an MTBDD for an n-qubit state |x>.
 * 
 * @param n Number of qubits.
 * @param x A qubitstring x \in {0,1}^n. 
 * 
 * @return An MTBDD encoding of the n-qubit state |x>.
 */
MTBDD mtbdd_create_basis_state_double(BDDVAR n, bool* x); //TODO: extend with complex_t
MTBDD mtbdd_create_basis_state_mpc(BDDVAR n, bool* x);

/**
 * Creates an MTBDD matrix which applies gate G to qubit t and I to all others.
 * 
 * @param n Total number of qubits.
 * @param t Target qubit for given gate (index connected first qubit to gate G).
 * @param I_dd mtbdd of predefined single qubit identity gate I.
 * @param G_dd mtbdd of predefined single qubit gate G.
 * 
 * @return An MTBDD encoding of I(0) x ... x G(t) x ... x I(n-1)
 * 
 */
MTBDD mtbdd_create_single_gate_for_qubits_mpc(BDDVAR n, BDDVAR t, MTBDD I_dd, MTBDD G_dd);

/**
 * Creates an MTBDD matrix which applies a control and gate G to qubit c and t and I to all others.
 * 
 * @param n Total number of qubits.
 * @param c Control qubit for given control gate
 * @param t Target qubit for given gate (index connected first qubit to gate G).
 * @param I_dd mtbdd of predefined single qubit identity gate I.
 * @param V00_dd mtbdd of |0><0| matrix
 * @param V11_dd mtbdd of |1><1| matrix
 * @param G_dd mtbdd of predefined single qubit gate G.
 * 
 * @return An MTBDD encoding of CU control unitary gate U={x,y,z, ....}
 * 
 */
MTBDD mtbdd_create_single_control_gate_for_qubits_mpc(BDDVAR n, BDDVAR c, BDDVAR t, MTBDD I_dd, MTBDD V00_dd, MTBDD V11_dd, MTBDD G_dd);

/**
 * Calculates the L2 norm of a mtbdd with leaves with mpc type.
 */
double mtbdd_getnorm_mpc(MTBDD dd, size_t nvars);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
