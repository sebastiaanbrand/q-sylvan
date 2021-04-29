#include "QASM_to_circuit.h"
#include "sylvan.h"

/**
 * Prints the final results of the measured qubits (flagged by <measurements>).
 * 
 * PARAMETERS:
 * - qdd: a QDD representing the statevector to be measured on
 * - measurements: an array containing booleans, where the to be measured qubits are flagged
 * - qubits: the number of qubits
 * - runs: the number of times the flagged qubits in the statevector should be measured
 * - show: toggle to print the measurement results
 */
void final_measure(QDD qdd, bool* measurements, BDDVAR qubits, BDDVAR shots, bool show);

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
 * Applies controlled <gate> to <qdd> on the qubit with index <i>.
 * 
 * PARAMETERS:
 * - qdd: the statevector on which to apply the controlled <gate>
 * - gate: the gate to be applied
 * - i: the index of the qubit in <qdd> on which to apply the controlled <gate>
 * 
 * RETURN:
 * - The statevector qdd where controlled <gate> has been applied on the qubit with index <i>
 */
#define apply_controlled_gate(qdd,gate,i,n) (CALL(apply_controlled_gate,qdd,gate,i,n));
TASK_DECL_4(QDD, apply_controlled_gate, QDD, Gate, BDDVAR, BDDVAR);

/**
 * Applies <gate> to <qdd> on the qubit with index <i>.
 * 
 * PARAMETERS:
 * - qdd: the statevector on which to apply <gate>
 * - gate: the gate to be applied
 * - i: the index of the qubit in <qdd> on which to apply <gate>
 * 
 * RETURN:
 * - The statevector qdd where <gate> has been applied on the qubit with index <i>
 */
#define apply_gate(qdd,gate,i,n) (CALL(apply_gate,qdd,gate,i,n));
TASK_DECL_4(QDD, apply_gate, QDD, Gate, BDDVAR, BDDVAR);

/**
 * Generates a gate QDD with <n> qubits that represents <gate> applied to qubit <k>
 * 
 * PARAMETERS:
 * - gate: the gate to be applied
 * - k: the index of the qubit on which to apply <gate>
 * - n: the total number of qubits
 * 
 * RETURN:
 * - The gate QDD with <n> qubits that represents <gate> applied to qubit <k>
 */
#define handle_control_matrix(gate,k,n) (CALL(handle_control_matrix,gate,k,n));
TASK_DECL_3(QDD, handle_control_matrix, Gate, BDDVAR, BDDVAR);

/**
 * Generates a gate QDD with <n> qubits that represents a swap operation between <qubit1> and <qubit2>
 * 
 * PARAMETERS:
 * - qubit1: the index of the first qubit in the swap
 * - qubit2: the index of the second qubit in the swap
 * - n: the total number of qubits
 * 
 * RETURN:
 * - The gate QDD with <n> qubits that represents a swap operation between <qubit1> and <qubit2>
 */
#define circuit_swap_matrix(qubit1,qubit2,n) (CALL(circuit_swap_matrix,qubit1,qubit2,n));
TASK_DECL_3(QDD, circuit_swap_matrix, BDDVAR, BDDVAR, BDDVAR);

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
#define check_classical_if(bits,expected,actual_list) (CALL(check_classical_if,bits,expected,actual_list));
TASK_DECL_3(bool, check_classical_if, BDDVAR, BDDVAR, int*);

/**
 * Runs the circuit_struct <c_s> for <runs> times and prints the results if <show> is true.
 * The run is done using the matrx-matrix method. Gate QDDs are multiplied with eachother
 * until the end of the circuit is reached, after which the resulting QDD is multiplied with
 * the statevector QDD. When number of nodes in the gate QDD exceed <limit>, the QDD is 
 * multiplied with the statevector QDD and the gate QDD is reset. This is to prevent exploding QDDs.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be run
 * - runs: the number of runs that should be done
 * - limit: the maximum number of nodes in the gate QDD
 * - show: prints the circuit results if set to true
 * 
 * RETURN:
 * - The resulting statevector QDD after running the circuit
 */
QDD run_c_struct_matrix(C_struct c_s, BDDVAR shots, bool show);

/**
 * Runs the circuit_struct <c_s> for <runs> times and prints the results if <show> is true.
 * The run is done using the matrx-vector method. Each gate is directly multiplied with
 * the statevector QDD.
 * 
 * PARAMETERS:
 * - c_s: the circuit_struct to be run
 * - runs: the number of runs that should be done
 * - show: prints the circuit results if set to true
 * 
 * RETURN:
 * - The resulting statevector QDD after running the circuit
 */
QDD run_c_struct(C_struct c_s, BDDVAR shots, bool show);