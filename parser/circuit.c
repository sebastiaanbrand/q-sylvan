#include "circuit.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

/*
 Create a circuit_struct with all attributes fitted to the circuit in the QASM file
 Attributes:
 - circuit: a two dimensional array containing gates corresponding to the QASM file
 - nvars: number of rows/qubits in the circuit
 - depth: number of columns in the circuit
 - qdd: A QDD state vector, set in the all-zero state
 - progress: a list of 'nvars' values that keep track the progress of each qubit, initialised to zero
 Input:
 - String: contains a path to a file describing a circuit in QASM
 Output:
 - circuit: struct containing all needed attributes
*/
Circuit* create_circuit(char* filename)
{
    C_struct c_s = make_c_struct(filename, true);
    Circuit *circuit_s = malloc(sizeof(Circuit));
    BDDVAR* progress = malloc(c_s.nvars * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s.nvars; i++) progress[i] = 0;
    circuit_s->progress = progress;
    circuit_s->nvars = c_s.nvars;
    circuit_s->depth = c_s.depth;
    circuit_s->qdd = qdd_create_all_zero_state(circuit_s->nvars);
    circuit_s->circuit = malloc(circuit_s->nvars*sizeof(circuit_s->circuit));
    for (BDDVAR i = 0; i < circuit_s->nvars; ++i) {
        circuit_s->circuit[i] = malloc(circuit_s->depth * sizeof(Gate));
        for (BDDVAR j = 0; j < circuit_s->depth; j++) {
            circuit_s->circuit[i][j] = c_s.circuit[i][j];
            circuit_s->circuit[i][j].control = malloc(c_s.circuit[i][j].controlSize*sizeof(BDDVAR));
            for (BDDVAR i = 0; i < c_s.circuit[i][j].controlSize; i++)
                circuit_s->circuit[i][j].control[i] = c_s.circuit[i][j].control[i];
        }
    }
    delete_c_struct(&c_s);
    return circuit_s;
}

/*
 Print the circuit represented by the circuit_struct
 Input:
 - circuit_s: struct that contains information about the circuit
 - vertical: the circuit is printed vertically if 'vertical' is true
 - show_rotation: show rotations for special rotation gates (e.g. Rx, Ry, Rz)
*/
void print_circuit(Circuit* circuit_s, bool vertical, bool show_rotation)
{
    bool has_rotation = false;
    bool negative_rotation = false;
    if (vertical) {
        for (BDDVAR j = 0; j < circuit_s->depth; j++) {
            for (BDDVAR i = 0; i < circuit_s->nvars; i++)
                printf(" |  ");
            printf("\n");
            for (BDDVAR i = 0; i < circuit_s->nvars; i++) {
                if (circuit_s->circuit[i][j].id == gate_barrier.id)
                    printf("----");
                else {
                    if (circuit_s->circuit[i][j].gateSymbol[0] == '-')
                        printf(" |");
                    else
                        printf(" %c",circuit_s->circuit[i][j].gateSymbol[0]);
                    if (circuit_s->circuit[i][j].gateSymbol[1] != '-')
                        printf("%c ",circuit_s->circuit[i][j].gateSymbol[1]);
                    else
                        printf("  ");
                }
            }
            printf("\n");
        }
    }
    else {
        for (BDDVAR i = 0; i < circuit_s->nvars; i++) {
            for (BDDVAR j = 0; j < circuit_s->depth; j++) {
                has_rotation = false;
                negative_rotation = false;
                for (BDDVAR k = 0; k < circuit_s->nvars; k++) {
                    if(circuit_s->circuit[k][j].id == 11 || circuit_s->circuit[k][j].id == 12 || circuit_s->circuit[k][j].id == 13) has_rotation = true;
                    if(circuit_s->circuit[k][j].rotation < 0) negative_rotation = true;
                }
                if (has_rotation && show_rotation) {
                    if(circuit_s->circuit[i][j].id == 11 || circuit_s->circuit[i][j].id == 12 || circuit_s->circuit[i][j].id == 13)
                        printf("-%s(%.4lf)",circuit_s->circuit[i][j].gateSymbol,roundf(circuit_s->circuit[i][j].rotation*10000)/10000);
                    else {
                        if (negative_rotation)
                            printf("-%s---------",circuit_s->circuit[i][j].gateSymbol);
                        else
                            printf("-%s--------",circuit_s->circuit[i][j].gateSymbol);
                    }
                }
                else
                    printf("-%s",circuit_s->circuit[i][j].gateSymbol);
            }
            printf("-\n");
        }
    }
}

/*
 Deletes the circuit_struct by freeing memory
 Input:
 - circuit_s: struct to be deleted
*/
void delete_circuit(Circuit* circuit_s)
{
    for (BDDVAR j = 0; j < circuit_s->depth; j++) {
        for (BDDVAR i = 0; i < circuit_s->nvars; i++) {
            if (circuit_s->circuit[i][j].id != gate_I.id) {
                if (circuit_s->circuit[i][j].controlSize > 0)
                    free(circuit_s->circuit[i][j].control);
            }
        }
    }
    for (BDDVAR i = 0; i < circuit_s->nvars; i++)
        free(circuit_s->circuit[i]);
    free(circuit_s->circuit);
    free(circuit_s->progress);
    free(circuit_s);
}

/*
 Apply in the qdd by applying a gate on the current progress step of qubit q and increment the progress counter
 Input:
 - circuit_s: struct containing information about the circuit to be advanced (advancement done in place)
 - q: qubit on which to advance a gate
 Output:
 - bool: true if gate is advanced, false if not possible to advance (e.g. a crtl is preceded by not yet advanced gates)
*/
TASK_IMPL_2(bool, advance, Circuit*, circuit_s, BDDVAR, q)
{
    Gate gate = circuit_s->circuit[q][circuit_s->progress[q]];
    BDDVAR gate_id = get_gateid(gate);
    for (BDDVAR i = 0; i < gate.controlSize; i++) {
        if (circuit_s->circuit[gate.control[i]][circuit_s->progress[gate.control[i]]].id != gate_ctrl.id)
            return false;
    }
    if (gate.controlSize == 0)
        circuit_s->qdd = qdd_gate(circuit_s->qdd, gate_id, q);
    else if (gate.controlSize == 1)
        circuit_s->qdd = qdd_cgate(circuit_s->qdd, gate_id, gate.control[0], q);
    else if (gate.controlSize == 2)
        circuit_s->qdd = qdd_cgate2(circuit_s->qdd, gate_id, gate.control[0], gate.control[1], q);
    else if (gate.controlSize == 3)
        circuit_s->qdd = qdd_cgate3(circuit_s->qdd, gate_id, gate.control[0], gate.control[1], gate.control[2], q);
    for (BDDVAR i = 0; i < gate.controlSize; i++)
        skip_gate(circuit_s, gate.control[i]);
    skip_gate(circuit_s, q);
    return true;
}

/*
 Helper function to get the qdd-based gateid from a gate_struct
 Input:
 - gate: struct from which to get the qdd-based gateid
 Output:
 - BDDVAR: the corresponding qdd-based gateid
*/
TASK_IMPL_1(BDDVAR, get_gateid, Gate, gate)
{
    BDDVAR gate_id;
    if (gate.id == gate_X.id)
        gate_id = GATEID_X;
    else if (gate.id == gate_Y.id)
        gate_id = GATEID_Y;
    else if (gate.id == gate_Z.id)
        gate_id = GATEID_Z;
    else if (gate.id == gate_H.id)
        gate_id = GATEID_H;
    else if (gate.id == gate_sX.id)
        gate_id = GATEID_sqrtX;
    else if (gate.id == gate_sY.id)
        gate_id = GATEID_sqrtY;
    else if (gate.id == gate_S.id)
        gate_id = GATEID_S;
    else if (gate.id == gate_Sd.id)
        gate_id = GATEID_Sdag;
    else if (gate.id == gate_T.id)
        gate_id = GATEID_T;
    else if (gate.id == gate_Td.id)
        gate_id = GATEID_Tdag;
    else if (gate.id == gate_Rx.id)
        gate_id = GATEID_Rx(gate.rotation);
    else if (gate.id == gate_Ry.id)
        gate_id = GATEID_Ry(gate.rotation);
    else if (gate.id == gate_Rz.id)
        gate_id = GATEID_Rz(gate.rotation);
    else {
        printf("Unknown gate: %d\n", gate.id);
        exit(1);
    }
    return gate_id;
}

/*
 Helper function to skip oncoming Identity and barrier gates after advancing a gate by advancing the progress counter
 Input:
 - circuit_s: circuit on which to apply the skips in the progress counter
 - i: the qubit on which to apply the skips in the progress counter
*/
void skip_gate(Circuit* circuit_s, BDDVAR i)
{
    bool depth, is_gate_I, is_gate_barrier;
    do {
        circuit_s->progress[i]++;
        depth = circuit_s->progress[i] < circuit_s->depth;
        is_gate_I = circuit_s->circuit[i][circuit_s->progress[i]].id == gate_I.id;
        is_gate_barrier = circuit_s->circuit[i][circuit_s->progress[i]].id == gate_barrier.id;
    } while (depth && (is_gate_I || is_gate_barrier));
}