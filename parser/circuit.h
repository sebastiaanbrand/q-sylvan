#include <stdbool.h>

#include "QASM_to_circuit.h"
#include "sylvan.h"

// Circuit struct
typedef struct Circuit {
    Gate** circuit;
    BDDVAR qubits;
    BDDVAR bits;
    BDDVAR depth;
    BDDVAR* progress;
    QDD qdd;
} Circuit;

/**
 * Create a circuit struct with all attributes fitted to the circuit in the QASM file.
 * 
 * ATTRIBUTES:
 * - circuit: a two dimensional array containing gates corresponding to the QASM file
 * - qubits: number of rows/qubits in the circuit
 * - depth: number of columns in the circuit
 * - qdd: A QDD state vector, set in the all-zero state
 * - progress: a list of 'qubits' values that keep track the progress of each qubit, initialised to zero
 * 
 * PARAMETERS:
 * - filename: contains a path to a file describing a circuit in QASM 
 * 
 * RETURN:
 * - circuit struct containing all needed attributes 
 */
Circuit* create_circuit(char* filename);

/**
 * Print the circuit represented by the circuit_struct
 * 
 * PARAMETERS:
 * - c_s: struct that contains information about the circuit
 * - vertical: the circuit is printed vertically if 'vertical' is true
 * - show_rotation: show rotations for special rotation gates (e.g. Rx, Ry, Rz)
 */
void print_circuit(Circuit* circuit_s, bool show_rotation);

/**
 * Deletes the circuit_struct by freeing memory
 * 
 * PARAMETERS:
 * - circuit_s: struct to be deleted
 */
void delete_circuit(Circuit* circuit_s);

/**
 * Skips oncoming identity and barrier gates after advancing a gate by advancing the progress counter
 * 
 * PARAMETERS:
 * - circuit_s: circuit on which to advance the progress counter
 * - i: the qubit on which to advance the progress counter
 * 
 * RETURN:
 * - true if successful, else false
 */
bool skip_to_gate(Circuit* circuit_s, BDDVAR i);

/**
 Apply in the qdd by applying a gate on the current progress step of qubit q and increment the progress counter

 * PARAMETERS:
 * - circuit_s: struct containing information about the circuit to be advanced
 * - q: qubit on which to advance a gate
 * 
 * RETURN:
 * - returns true if gate is advanced, false if not possible to advance (e.g. a crtl is preceded by not yet advanced gates)
 */
#define advance(circuit_s,q) (CALL(advance,circuit_s,q));
TASK_DECL_2(bool, advance, Circuit*, BDDVAR);

/**
 Helper function to get the qdd-based gateid from a gate_struct

 * PARAMETERS:
 * - gate: struct from which to get the qdd-based gateid
 * 
 * RETURN:
 * - returns the corresponding qdd-based gateid
 */
#define get_gateid(gate) (CALL(get_gateid,gate));
TASK_DECL_1(BDDVAR, get_gateid, Gate);
