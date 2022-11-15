#include "QASM_to_circuit.h"

/**
 * Prints the final results of the measured qubits (flagged by <measurements>).
 * 
 * PARAMETERS:
 * - qmdd: a QMDD representing the statevector to be measured on
 * - measurements: an array containing booleans, where the to be measured qubits are flagged
 * - qubits: the number of qubits
 * - runs: the number of times the flagged qubits in the statevector should be measured
 * - show: toggle to print the measurement results
 */
void final_measure(QMDD qmdd, int* measurements, C_struct c_s, bool* results);

/**
 * Returns the Sylvan GATEID that corresponds to <gate>
 * 
 * PARAMETERS:
 * - gate: gate_struct of which to get the Sylvan GATEID
 * 
 * RETURN:
 * - The Sylvan GATEID which corresponds to <gate>
 */
#define get_gate_id(gate) (CALL(get_gate_id,gate));
TASK_DECL_1(BDDVAR, get_gate_id, Gate);

/**
 * Applies controlled <gate> to <qmdd> on the qubit with index <i>.
 * 
 * PARAMETERS:
 * - qmdd: the statevector on which to apply the controlled <gate>
 * - gate: the gate to be applied
 * - i: the index of the qubit in <qmdd> on which to apply the controlled <gate>
 * 
 * RETURN:
 * - The statevector qmdd where controlled <gate> has been applied on the qubit with index <i>
 */
#define apply_controlled_gate(qmdd,gate,i,n) (CALL(apply_controlled_gate,qmdd,gate,i,n));
TASK_DECL_4(QMDD, apply_controlled_gate, QMDD, Gate, BDDVAR, BDDVAR);

/**
 * Applies <gate> to <qmdd> on the qubit with index <i>.
 * 
 * PARAMETERS:
 * - qmdd: the statevector on which to apply <gate>
 * - gate: the gate to be applied
 * - i: the index of the qubit in <qmdd> on which to apply <gate>
 * 
 * RETURN:
 * - The statevector qmdd where <gate> has been applied on the qubit with index <i>
 */
#define apply_gate(qmdd,gate,i,n) (CALL(apply_gate,qmdd,gate,i,n));
TASK_DECL_4(QMDD, apply_gate, QMDD, Gate, BDDVAR, BDDVAR);

/**
 * Generates a gate QMDD with <n> qubits that represents <gate> applied to qubit <k>
 * 
 * PARAMETERS:
 * - gate: the gate to be applied
 * - k: the index of the qubit on which to apply <gate>
 * - n: the total number of qubits
 * 
 * RETURN:
 * - The gate QMDD with <n> qubits that represents <gate> applied to qubit <k>
 */
#define handle_control_matrix(gate,k,n) (CALL(handle_control_matrix,gate,k,n));
TASK_DECL_3(QMDD, handle_control_matrix, Gate, BDDVAR, BDDVAR);

/**
 * Generates a gate QMDD with <n> qubits that represents a swap operation between <qubit1> and <qubit2>
 * 
 * PARAMETERS:
 * - qubit1: the index of the first qubit in the swap
 * - qubit2: the index of the second qubit in the swap
 * - n: the total number of qubits
 * 
 * RETURN:
 * - The gate QMDD with <n> qubits that represents a swap operation between <qubit1> and <qubit2>
 */
#define circuit_swap_matrix(qubit1,qubit2,n) (CALL(circuit_swap_matrix,qubit1,qubit2,n));
TASK_DECL_3(QMDD, circuit_swap_matrix, BDDVAR, BDDVAR, BDDVAR);

/**
 * Checks if the measuring gate at position (<qubit>,<depth>) in <c_s> is the last gate on that qubit
 * 
 * PARAMETERS:
 * - c_s: the circuit on which to check this 
 * - qubit: the index of the qubit
 * - depth: the depth of the gate
 * 
 * RETURN:
 * - true if it is the last gate, else false
 */
#define is_final_measure(c_s,qubit,depth) (CALL(is_final_measure,c_s,qubit,depth));
TASK_DECL_3(bool, is_final_measure, C_struct, BDDVAR, BDDVAR);

/**
 * Checks if the classical register contains the expected value of the classical if statement
 * 
 * PARAMETERS:
 * - bits: the number of bits in the register
 * - expected: the expected number to see in the register
 * - actual_list: the register containing the bits
 * 
 * RETURN:
 * - true if they are equal, else false
 */
#define check_classical_if(bits,gate,actual_list) (CALL(check_classical_if,bits,gate,actual_list));
TASK_DECL_3(bool, check_classical_if, BDDVAR, Gate, bool*);

/**
 * Checks if the circuit has intermediate measuring gates, or if the circuit is only
 * measured at the end of the circuit.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be checked
 * 
 * RETURN:
 * - true if there are intermediate  measuring gates, else false
 */
bool check_measuring_gates(C_struct c_s);

/**
 * Advances the progress of qubit <i> to the next gate
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct
 * - progress: the progress array
 * - i: the qubit index
 * - skip_barrier: skips barrier if true
 */
void skip(C_struct c_s, BDDVAR* progress, BDDVAR i, bool skip_barrier);

/**
 * Measures qubit <i>.
 * 
 * PARAMETERS:
 * - qmdd: the qmdd to be measured on
 * - i: the qubit to measure
 * - n: the number of qubits in <qmdd>
 * - result: the result of the measurement
 * 
 * RETURN:
 * - the resulting qmdd after measuring
 **/
QMDD measure(QMDD qmdd, BDDVAR i, BDDVAR n, bool* result);

/**
 * Checks for the greedy method if the controls of a gate are free and if it satisfies the classical
 * control (if it is classically controlled).
 * 
 * PARAMETERS:
 * - c_s: the circuit struct
 * - gate: the gate to check
 * - progress: the progress counter of all qubits
 * - i: the gate index
 * - results: the classical register
 * 
 * RETURN:
 * - true if everything is satisfied, else false
 **/
bool check_gate(C_struct c_s, Gate gate, BDDVAR* progress, BDDVAR i, bool* results);

/**
 * A subfunction which runs the greedy method from <column> to the next barrier.
 * 
 * PARAMETERS:
 * - c_s: the circuit struct
 * - prev_qmdd: the current qmdd
 * - column: the current column
 * - measurements: the final measurement flags
 * - results: the classical register
 * - experiments: prints nodecounts if true
 * - n_gates: the current number of nodes already applied
 * 
 * RETURN:
 * - the resulting qmdd after applying all gates between column and the next barrier
 **/
QMDD greedy(C_struct c_s, QMDD prev_qmdd, BDDVAR* column, int* measurements, bool* results, bool experiments, BDDVAR* n_gates);

/**
 * A subfunction which runs the matmat method from <column> to the next barrier.
 * 
 * PARAMETERS:
 * - c_s: the circuit struct
 * - vec: the current qmdd
 * - column: the current column
 * - measurements: the final measurement flags
 * - results: the classical register
 * - limit: the node-threshold to check if a matvec multiplication should be done
 * - experiments: prints nodecounts if true
 * - n_gates: the current number of nodes already applied
 * 
 * RETURN:
 * - the resulting qmdd after applying all gates between column and the next barrier
 **/
QMDD matmat(C_struct c_s, QMDD vec, BDDVAR* column, int* measurements, bool* results, int limit, bool experiments, BDDVAR* n_gates);

/**
 * Runs the circuit_struct <c_s> with a greedy method and stores measurements in <bit_res>. 
 * The run is done using the matrx-vector method. Each gate is directly multiplied with 
 * the statevector QMDD.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be run
 * - bit_res: the bit register to store measurements in
 * - nodes: the maximum number of nodes in the QMDD at any point in the algorithm
 * 
 * RETURN:
 * - The resulting statevector QMDD after running the circuit
 */
QMDD greedy_run_circuit(C_struct c_s, int* measurements, bool* results, bool experiments);

/**
 * Runs the circuit_struct <c_s> for <runs>  and stores measurements in <bit_res>.
 * The run is done using the matrix-matrix method. Gate QMDDs are multiplied with eachother
 * until the end of the circuit is reached, after which the resulting QMDD is multiplied with
 * the statevector QMDD. When number of nodes in the gate QMDD exceed <limit>, the QMDD is 
 * multiplied with the statevector QMDD and the gate QMDD is reset. This is to prevent exploding QMDDs.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be run
 * - bit_res: the bit register to store measurements in
 * - limit: the maximum number of nodes in the gate QMDD
 * - nodes: the maximum number of nodes in the QMDD at any point in the algorithm
 * 
 * RETURN:
 * - The resulting statevector QMDD after running the circuit
 */
QMDD run_c_struct_matrix(C_struct c_s, int* measurements, bool* bit_res, int limit, bool experiments);

/**
 * Runs the circuit_struct <c_s> and stores measurements in <bit_res>. The run is done using the
 * matrix-vector method. Each gate is directly multiplied with the statevector QMDD.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be run
 * - bit_res: the bit register to store measurements in
 * - nodes: the maximum number of nodes in the QMDD at any point in the algorithm
 * 
 * RETURN:
 * - The resulting statevector QMDD after running the circuit
 */
QMDD run_c_struct(C_struct c_s, int* measurements, bool* bit_res, bool experiments);