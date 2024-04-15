#ifndef QSYLVAN_SIMULATOR_MTBDD_H
#define QSYLVAN_SIMULATOR_MTBDD_H

#include <sylvan_int.h>
#include <qsylvan_gates.h>


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

static const uint64_t MTBDD_TERMINAL = 1;
static const BDDVAR   MTBDD_INVALID_VAR = UINT8_MAX;

/**
 * Creates an MTBDD for an n-qubit state |00...0>.
 * 
 * @param n Number of qubits.
 * 
 * @return An MTBDD encoding the n-qubit state |00..0>.
 */
MTBDD mtbdd_create_all_zero_state(BDDVAR n);

/**
 * Creates an MTBDD for an n-qubit state |x>.
 * 
 * @param n Number of qubits.
 * @param x A bitstring x \in {0,1}^n. 
 * 
 * @return An MTBDD encoding of the n-qubit state |x>.
 */
MTBDD mtbdd_create_basis_state(BDDVAR n, bool* x);

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
*/
MTBDD
mtbdd_stack_matrix(BDDVAR k, gate_id_t gateid);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
