#ifdef __cplusplus
 extern "C"
{
#endif

typedef enum operation_type {
    op_gate,
    op_measurement,
    op_blank
} operation_type_t;

typedef struct quantum_op_s {
    operation_type_t type;
    char name[32];
    double angle[3];
    int targets[2];
    int ctrls[3];
    int meas_dest;
    struct quantum_op_s* next;
} quantum_op_t;

typedef struct quantum_circuit_s {
    int qreg_size;
    int creg_size;
    bool *creg;
    bool has_intermediate_measurements;
    quantum_op_t *operations;
} quantum_circuit_t;

/**
 * quantum_op_t is a linked-list where every element is malloc'ed.
 * free_quantum_ops() should be called on this after use.
 */
quantum_circuit_t* parse_qasm_file(char *filepath);


/**
 * Invert qubit order if that yields less controls below target qubits.
*/
void optimize_qubit_order(quantum_circuit_t *circuit);

/**
 * free()'s all quantum_ops from and including 'first'.
*/
void free_quantum_circuit(quantum_circuit_t* circuit);

void print_quantum_op(quantum_op_t *op);
void print_quantum_circuit(quantum_circuit_t* circuit);
void print_creg(quantum_circuit_t* circuit);

#ifdef __cplusplus
}
#endif