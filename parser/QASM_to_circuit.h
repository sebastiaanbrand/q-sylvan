#include <stdbool.h>

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

struct gate {
    BDDVAR id;
    char gateSymbol[2];
    float rotation;
    BDDVAR gateId;
} gate;

static const struct gate gate_I = {0, "--", 0, GATEID_I};
static const struct gate gate_X = {1, "X-", 0.5, GATEID_X};
static const struct gate gate_Y = {2, "Y-", 0.5, GATEID_Y};
static const struct gate gate_Z = {3, "Z-", 0.5, GATEID_Z};
static const struct gate gate_H = {4, "H-", 0, GATEID_H};
static const struct gate gate_S = {5, "S-", 0.25, GATEID_S};
static const struct gate gate_Sd = {6, "Sd", -0.25, GATEID_Sdag};
static const struct gate gate_T = {7, "T-", 0.125, GATEID_T};
static const struct gate gate_Td = {8, "Td", -0.125, GATEID_Tdag};
static const struct gate gate_sX = {9, "sX", 0.25, GATEID_sqrtX};
static const struct gate gate_sY = {10, "sY", 0.25, GATEID_sqrtY};
static const struct gate gate_Rx = {11, "Rx", 0, 0};
static const struct gate gate_Ry = {12, "Ry", 0, 0};
static const struct gate gate_Rz = {13, "Rz", 0, 0};
static const struct gate barrier = {14, "|-", 0, 0};

void make_circuit(char *filename, bool optimize);
void print_circuit(int32_t** circuit, BDDVAR* nvars, BDDVAR* curr_depth);
BDDVAR optimize_circuit(int32_t** circuit, BDDVAR nvars, BDDVAR curr_depth);
void optimize_controlled_gate(int32_t** circuit, BDDVAR depth1, BDDVAR depth2, BDDVAR target, BDDVAR nvars);
void find_palindromes(int32_t** circuit, BDDVAR curr_depth, BDDVAR var, BDDVAR nvars, BDDVAR depth);
BDDVAR reduce_circuit(int32_t** circuit, BDDVAR nvars, BDDVAR depth);
void reduce_controlled_gate(int32_t** circuit, BDDVAR depth, BDDVAR target, BDDVAR nvars);

#define get_qubits_circuit(token,n_qubits,qubits) (CALL(get_qubits_circuit,token,n_qubits,qubits));
TASK_DECL_3(bool, get_qubits_circuit, char*, BDDVAR, BDDVAR*);

#define get_gateid_circuit(tokens,gate_id) (CALL(get_gateid_circuit,tokens,gate_id));
TASK_DECL_2(bool, get_gateid_circuit, char*, uint32_t*);

#define get_parallel_depth(targets,circuit,n_qubits,curr_depth,gateid) (CALL(get_parallel_depth,targets,circuit,n_qubits,curr_depth,gateid));
TASK_DECL_5(bool, get_parallel_depth, char*, int32_t**, BDDVAR, BDDVAR*, uint32_t*);

#define handle_line_circuit(line,circuit,nvars,curr_depth) (CALL(handle_line_circuit,line,circuit,nvars,curr_depth));
TASK_DECL_4(bool, handle_line_circuit, char*, int32_t**, BDDVAR*, BDDVAR*);