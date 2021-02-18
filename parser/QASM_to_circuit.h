#include <stdbool.h>

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

typedef struct gate Gate;
struct gate
{
    BDDVAR id;
    char gateSymbol[2];
    float rotation;
    BDDVAR *control;
    BDDVAR controlSize;
};

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

static const BDDVAR max_qubits = 128;
static const BDDVAR max_wire = 1024;
static BDDVAR wire_i;
static BDDVAR* nvars;
static BDDVAR* curr_depth;

void circuit_exit(Gate** circuit);
Gate** make_circuit(char *filename, bool optimize);
void print_circuit(Gate** circuit, BDDVAR* nvars, BDDVAR* curr_depth);
BDDVAR optimize_circuit(Gate** circuit, BDDVAR nvars, BDDVAR curr_depth);
void optimize_circuit_p(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2, BDDVAR nvars, BDDVAR depth);
void optimize_controlled_gate(Gate** circuit, BDDVAR depth1, BDDVAR depth2, BDDVAR target, BDDVAR nvars);
bool find_palindromes(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2);
void remove_gates(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2);
BDDVAR reduce_circuit(Gate** circuit, BDDVAR nvars, BDDVAR depth);
void reduce_gate(Gate** circuit, BDDVAR target, BDDVAR depth);
BDDVAR get_reduce_depth(Gate** circuit, BDDVAR target, BDDVAR depth);

bool get_qubits_circuit(char* token, BDDVAR n_qubits, BDDVAR* qubits);
void copy_Gate(Gate src, Gate* dst);
BDDVAR get_gateid_circuit(char* gate_str, Gate* gate);
bool handle_gate(char* targets, Gate** circuit, BDDVAR n_qubits, BDDVAR* curr_depth, Gate Gate);
bool handle_line_circuit(char* line, Gate** circuit, BDDVAR* nvars, BDDVAR* curr_depth);
