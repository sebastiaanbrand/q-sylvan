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

#include <sylvan_int.h>
#include <qsylvan_gates.h>


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
MTBDD mtbdd_create_all_zero_state_double(BDDVAR n);
MTBDD mtbdd_create_all_zero_state_mpc(BDDVAR n);

/**
 * Creates an MTBDD for an n-qubit state |x>.
 * 
 * @param n Number of qubits.
 * @param x A qubitstring x \in {0,1}^n. 
 * 
 * @return An MTBDD encoding of the n-qubit state |x>.
 */
MTBDD mtbdd_create_basis_state_double(BDDVAR n, bool* x);
MTBDD mtbdd_create_basis_state_mpc(BDDVAR n, bool* x);

/**
 * Creates an MTBDD matrix which applies gate U to qubit t and I to all others.
 * 
 * @param n Total number of qubits.
 * @param t Target qubit for given gate.
 * @param gateid Gate ID of predefined single qubit gate U.
 * 
 * @return An MTBDD encoding of I_0 \tensor ... U_t ... \tensor I_{n-1}
 */
MTBDD mtbdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid);

/**
 *
 * 
 * 
 */
MTBDD
mtbdd_stack_matrix(BDDVAR k, gate_id_t gateid);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
