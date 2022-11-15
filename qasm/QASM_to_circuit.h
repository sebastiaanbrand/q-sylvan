#include <stdbool.h>

#include "qsylvan.h"

// Gates struct
typedef struct Gate
{
    BDDVAR id;
    char gateSymbol[2];
    float rotation;
    BDDVAR *control;
    BDDVAR controlSize;
    int16_t classical_control;
    int16_t classical_expect;
} Gate;

// Default gates
static const Gate gate_I = {0, "--", 0, NULL, 0, -1, -1};
static const Gate gate_X = {1, "X-", 0.5, NULL, 0, -1, -1};
static const Gate gate_Y = {2, "Y-", 0.5, NULL, 0, -1, -1};
static const Gate gate_Z = {3, "Z-", 0.5, NULL, 0, -1, -1};
static const Gate gate_H = {4, "H-", 0.5, NULL, 0, -1, -1};
static const Gate gate_sX = {9, "sX", 0.25, NULL, 0, -1, -1};
static const Gate gate_sY = {10, "sY", 0.25, NULL, 0, -1, -1};
static const Gate gate_S = {5, "S-", 0.25, NULL, 0, -1, -1};
static const Gate gate_Sd = {6, "Sd", -0.25, NULL, 0, -1, -1};
static const Gate gate_T = {7, "T-", 0.125, NULL, 0, -1, -1};
static const Gate gate_Td = {8, "Td", -0.125, NULL, 0, -1, -1};
static const Gate gate_Rx = {11, "Rx", 0, NULL, 0, -1, -1};
static const Gate gate_Ry = {12, "Ry", 0, NULL, 0, -1, -1};
static const Gate gate_Rz = {13, "Rz", 0, NULL, 0, -1, -1};
static const Gate gate_ctrl = {14, "@-", 0, NULL, 0, -1, -1};
static const Gate gate_ctrl_c = {15, "|-", 0, NULL, 0, -1, -1};
static const Gate gate_measure = {16, "M-", 0, NULL, 0, -1, -1};
static const Gate gate_barrier = {17, "#-", 0, NULL, 0, -1, -1};

// Circuit struct
typedef struct C_struct
{
    Gate** circuit;
    BDDVAR qubits;
    BDDVAR bits;
    BDDVAR depth;
    BDDVAR max_qubits;
    BDDVAR max_wire;
} C_struct;

// Default circuit
static const C_struct c_struct_default = {NULL, 0, 0, 0, 128, 1024};

/**
 * Creates a circuit struct that represents the circuit described in <filename>. If
 * <optimize> is true, negating gates will be removed.
 * 
 * ATTRIBUTES:
 * - circuit: a two dimensional array containing gates corresponding to the QASM file
 * - qubits: number of rows/qubits in the circuit
 * - depth: number of columns in the circuit
 * - max_qubits: the maximum number of qubits
 * - max_wire: the maximum wire length (can be extended with 'reallocate_wire')
 * 
 * PARAMETERS:
 * - filename: the path to the file containing the QASM circuit code
 * - optimize: if optimize is true, negating gates will be removed
 * 
 * RETURN:
 * - a circuit struct reperesenting the circuit described in <filename>
 * 
 * NOTE:
 * The current maximum number of qubits is 128. However this can easily be changed.
 */
C_struct make_c_struct(char *filename, bool optimize);

/**
 * Increases the maximum depth of <c_s> and reallocates the circuit.
 * 
 * PARAMETERS:
 * - c_s: circuit struct of which the depth needs to be increased
 */
void reallocate_wire(C_struct* c_s);

/**
 * Returns a copy of <c_s> where max_qubits and max_wire are set to qubits and depth respectively.
 * 
 * PARAMETERS:
 * - c_s: circuit struct to copy
 * 
 * RETURN:
 * - A new circuit struct which is a copy of the original (but max_qubits and max_wire changed)
 */
C_struct copy_c_struct(C_struct* c_s);

/**
 * Frees all variables used by <c_s>
 * 
 * PARAMETERS:
 * - c_s: the circuit struct to free
 */
void delete_c_struct(C_struct* c_s);

/**
 * Parses a gate command <line> to be included in <c_s>
 * 
 * PARAMETERS:
 * - line: the gate command to be included
 * - c_s: the circuit struct in which the gate should be placed
 * 
 * RETURN:
 * - true if succeeded, else false
 */
bool handle_line_c_struct(char* line, C_struct* c_s);

/**
 * Fills <qubits> with the indices of <n_qubits> qubits described in <token>
 * 
 * PARAMETERS:
 * - token: the string that contains the indices of all qubits
 * - n_qubits: the number of indices to find in <token>
 * - qubits: the array where the indices should be stored
 * 
 * RETURN:
 * - true if succeeded, false otherwise.
 */
bool get_qubits_c_struct(char* token, BDDVAR n_qubits, BDDVAR* qubits);

/**
 * Get the corresponding gate based on <gate_str>. Also count how many controls the gate has (including target).
 * 
 * PARAMETERS:
 * - gate_str: the string which contains QASM code describing a gate
 * - gate_s: a pointer to a gate location where the resulting gate should be stored
 * 
 * RETURN:
 * - the number of controls (including target) <gate_s> has. Aborts if the gate is unknown.
 */
int get_gateid_c_struct(char* gate_str, Gate* gate);

/**
 * Stores a measurement gate.
 * 
 * PARAMETERS:
 * - c_s: the circuit in which to store the measurement gate
 * - gate_s: the measurement gate struct (needed so we do not overwrite the default?)
 * - qubits: the target qubit for the gate and target bit for the measurement
 */
void handle_measure(C_struct* c_s, Gate gate_s, BDDVAR* qubits);

/**
 * Stores a barrier gate across the whole current column.
 * 
 * PARAMETERS:
 * - c_s: the circuit in which to store the barrier column
 * - gate_s: the barrier gate struct (needed so we do not overwrite the default?)
 */
void handle_barrier(C_struct* c_s, Gate gate_s);

/**
 * Swap two array elements
 * 
 * PARAMETERS:
 * - a: the first array element
 * - b: the second array element
 */
void swap(BDDVAR* a, BDDVAR* b);

/**
 * A partition function for quicksort
 * 
 * PARAMETERS:
 * - qubits: the circuit in which to store the barrier column
 * - low: the first index of the array
 * - high: the last index of  the array
 * 
 * RETURN:
 * - the pivot index
 */
BDDVAR partition(BDDVAR* qubits, int low, int high);

/**
 * Uses quicksort to sort the qubits
 * 
 * PARAMETERS:
 * - qubits: the circuit in which to store the barrier column
 * - low: the first index of the array
 * - high: the last index of  the array
 */
void sort(BDDVAR* qubits, int low, int high);

/**
 * Stores <gate_s> and corresponding control gates on all indices of <qubits> in <c_s>.
 * 
 * PARAMETERS:
 * - c_s: the circuit where everything should be stored
 * - n_qubits: the size of <qubits> - 1 (minus the target index)
 * - qubits: an array containing the indices of all control qubits (and target qubit)
 * - gate_s: the gate struct to be stored in <c_s>
 */
void handle_gate(C_struct* c_s, BDDVAR n_qubits, BDDVAR* qubits, Gate Gate);

/**
 * Gets the depth of the first preceding (or successive) gate of the gate at position (<q>,<depth>) in 
 * <c_s> if <successive> is false (or true).
 * 
 * PARAMETERS:
 * - c_s: the circuit to be used
 * - q: the qubit the gate is on
 * - depth: the column of the current gate
 * - successive: a variable which tells if the gate should be preceding or successive
 * 
 * RETURN:
 * - the depth of the first preceding (or successive) gate
 */
BDDVAR get_next_gate(C_struct* c_s, BDDVAR q, BDDVAR depth, bool successive);

/**
 * Optimize <c_s> by removing redundant gates by finding a set of gates which is an even palindrome 
 * of negating gates, e.g. --Z--RX(0.1)--H--H--RX(-0.1)--Z--
 * 
 * PARAMETERS:
 * - c_s: the circuit struct to optimize
 */
void optimize_c_struct(C_struct* c_s);

/**
 * Helper function for optimizing to be called recursively (since we dont know palindrome length) remove negating gates
 * 
 * PARAMETERS:
 * - c_s: the circuit struct to optimize
 * - q: the qubit we are looking at while finding palindromes
 * - depth1: the column containing the first gate of the negating pair
 * - depth2: the column containing the second gate of the negating pair
 */
void optimize_c_struct_p(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2);

/**
 * Checks in <c_s> if the gate on qubit <q> at <depth1> is the same as <depth2>.
 * 
 * PARAMETERS:
 * - c_s: the circuit struct to optimize
 * - q: the qubit we are looking at
 * - depth1: the column containing the first gate
 * - depth2: the column containing the second gate
 * 
 * RETURN:
 * - returns true if the two gates are a negation and there is nothing in between and false otherwise
 */
bool find_palindromes(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2);

/**
 * Removes a gates from <c_s> on qubit <q> at <depth>.
 * 
 * PARAMETERS:
 * - c_s: the circuit where the gate should be removed
 * - q: the qubit on which the gate is
 * - depth: the column of the gate
 */
void remove_gate(C_struct* c_s, BDDVAR q, BDDVAR depth);

/**
 * Reduces the depth variable of <c_s>. First all gates are moved to the left if possible. (parallelising)
 * Then the depth variable is set to the rightmost gate that is not identity.
 * 
 * PARAMETERS:
 * - c_s: the circuit which depth should be reduced
 */
void reduce_c_struct(C_struct* c_s);

/**
 * Reduces the depth of a gate on qubit <target> in column <depth> which is classically controlled.
 * 
 * PARAMETERS:
 * - c_s: the circuit in which the depth a gate needs to be reduced
 * - target: the qubit on which the gate is
 * - depth: the column in which the gate is
 */
void reduce_classical_gate(C_struct* c_s, BDDVAR target, BDDVAR depth);

/**
 * Reduces the depth of a gate on qubit <target> in column <depth>.
 * 
 * PARAMETERS:
 * - c_s: the circuit in which the depth a gate needs to be reduced
 * - target: the qubit on which the gate is
 * - depth: the column in which the gate is
 */
void reduce_gate(C_struct* c_s, BDDVAR target, BDDVAR depth);

/**
 * Given the gate in <c_s> on qubit <target> at column <depth>, get the depth it can be reduced to.
 * 
 * PARAMETERS:
 * - c_s: the circuit containing the gate
 * - target: the qubit on which the gate lies
 * - depth: the column where the gate lies
 * 
 * RETURN:
 * - the reducable depth for the gate
 */
BDDVAR get_reduce_depth(C_struct* c_s, BDDVAR target, BDDVAR depth);

/**
 * Prints the circuit contained in <c_s>. If <show_rotation> the special rotation gates will have their rotations printed too.
 * 
 * PARAMETERS:
 * - c_s: the circuit struct that contains circuit to be printed
 * - show_rotation: if true, the special rotation gates will have their rotations printed
 */
void print_c_struct(C_struct c_s, bool show_rotation);
