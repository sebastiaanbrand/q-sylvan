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

#ifdef __cplusplus
 extern "C"
{
#endif


/**
 * The structure quantum_circuit_t contains the QASM quantum circuit as defined in the QASM file.
 * 
 * It is composed of the structure quantum_op_t which is a linked-list of gates or gate-control 
 * quantum subcircuits of the quantum circuit.
 * 
 * Every element is put on the heap (malloc'ed) and should be released after use with free_quantum_ops().
 *
 */
typedef enum operation_type {

    op_gate,                            // gate element
    op_measurement,                     // measurement element
    op_blank                            // empty element

} operation_type_t;

typedef struct quantum_op_s {

    operation_type_t type;               // Element can be a gate or a measurement or blank
    char name[32];                       // QASM name of this gate
    double angle[3];                     // Arguments for Rotation gates(theta) or U(theta, phi, lambda), could be empty = 0,0,0, double -> mpfr, after library parser
    int targets[2];                      // Index of qubit which is connected to the first or second non-identity gate 
    int ctrls[3];                        // Index of qubit which is connected to the first, second or third control o
    int meas_dest;                       // Measurement dest
    struct quantum_op_s* next;

} quantum_op_t;

typedef struct quantum_circuit_s {

    char name[200];                       // Circuit name
    int qreg_size;                        // qubit register size
    int creg_size;                        // classical register size
    bool *creg;                           // bool is used as bit here, classical register ( size = sizeof(bool)*creg_size )
    bool has_intermediate_measurements;   // If true measurements not only at the end
    quantum_op_t *operations;             // First element of linked-list with all gates or gate-coontrol pairs
    bool reversed_qubit_order;            // Order of indices of the qubits, if true numbering is from bottom to top

} quantum_circuit_t;

/**
 * Parser of the QASM file, returns the quantum circuit in above structures.
 */
quantum_circuit_t* parse_qasm_file(char *filepath);

/**
 * Invert qubit order if that yields less controls below target qubits.
 */
void optimize_qubit_order(quantum_circuit_t *circuit);

/**
 * Free all quantum elements found in quantum_op_s including 'first'. What is 'first'? 
 */
void free_quantum_circuit(quantum_circuit_t* circuit);

/**
 * Print operations of quantum circuit, quantum elements or control registers.
 */
void print_quantum_op(quantum_op_t *op);
void print_quantum_circuit(quantum_circuit_t* circuit);
void fprint_creg(FILE *stream, quantum_circuit_t* circuit);


#ifdef __cplusplus
}
#endif