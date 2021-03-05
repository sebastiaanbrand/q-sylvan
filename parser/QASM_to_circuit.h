#include <stdbool.h>

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

// Gates struct
typedef struct Gate
{
    BDDVAR id;
    char gateSymbol[2];
    float rotation;
    BDDVAR *control;
    BDDVAR controlSize;
} Gate;

// Default gates
static const Gate gate_I = {0, "--", 0, NULL, 0};
static const Gate gate_X = {1, "X-", 0.5, NULL, 0};
static const Gate gate_Y = {2, "Y-", 0.5, NULL, 0};
static const Gate gate_Z = {3, "Z-", 0.5, NULL, 0};
static const Gate gate_H = {4, "H-", 0.5, NULL, 0};
static const Gate gate_sX = {9, "sX", 0.25, NULL, 0};
static const Gate gate_sY = {10, "sY", 0.25, NULL, 0};
static const Gate gate_S = {5, "S-", 0.25, NULL, 0};
static const Gate gate_Sd = {6, "Sd", -0.25, NULL, 0};
static const Gate gate_T = {7, "T-", 0.125, NULL, 0};
static const Gate gate_Td = {8, "Td", -0.125, NULL, 0};
static const Gate gate_Rx = {11, "Rx", 0, NULL, 0};
static const Gate gate_Ry = {12, "Ry", 0, NULL, 0};
static const Gate gate_Rz = {13, "Rz", 0, NULL, 0};
static const Gate gate_ctrl = {14, "@-", 0, NULL, 0};
static const Gate gate_measure = {15, "M-", 0, NULL, 0};
static const Gate gate_barrier = {16, "|-", 0, NULL, 0};

// Circuit struct
typedef struct C_struct
{
    Gate** circuit;
    BDDVAR nvars;
    BDDVAR depth;
    BDDVAR max_qubits;
    BDDVAR max_wire;
} C_struct;

// Default circuit
static const C_struct c_struct_default = {NULL, 0, 0, 128, 1024};

// Function definitions
C_struct make_c_struct(char *filename, bool optimize);
void reallocate_wire(C_struct* c_s);
void delete_c_struct(C_struct* c_s);
bool handle_line_c_struct(char* line, C_struct* c_s);
bool get_qubits_c_struct(char* token, BDDVAR n_qubits, BDDVAR* qubits);
BDDVAR get_gateid_c_struct(char* gate_str, Gate* gate);
void copy_Gate(Gate src, Gate* dst);
void handle_barrier(C_struct* c_s, Gate gate_s);
void handle_gate(C_struct* c_s, BDDVAR n_qubits, BDDVAR* qubits, Gate Gate);
void optimize_c_struct(C_struct* c_s);
void optimize_c_struct_p(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2);
bool find_palindromes(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2);
void remove_gates(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2);
void reduce_c_struct(C_struct* c_s);
void reduce_gate(C_struct* c_s, BDDVAR target, BDDVAR depth);
BDDVAR get_reduce_depth(C_struct* c_s, BDDVAR target, BDDVAR depth);
void print_c_struct(C_struct c_s, bool vertical, bool show_rotation);
