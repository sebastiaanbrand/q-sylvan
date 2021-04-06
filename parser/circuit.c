#include "circuit.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

Circuit* create_circuit(char* filename)
{
    // Create a c struct based on the file
    C_struct c_s = make_c_struct(filename, true);
    // Create a circuit struct
    Circuit *circuit_s = malloc(sizeof(Circuit));
    // Create a progress array
    BDDVAR* progress = malloc(c_s.nvars * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s.nvars; i++) progress[i] = 0;

    // Copy important info
    circuit_s->progress = progress;
    circuit_s->nvars = c_s.nvars;
    circuit_s->depth = c_s.depth;
    // Create an initial all-zero qdd
    circuit_s->qdd = qdd_create_all_zero_state(circuit_s->nvars);
    // Copy the circuit
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
    // Remove the c struct (since everything we need is copied)
    delete_c_struct(&c_s);
    return circuit_s;
}

void print_circuit(Circuit* circuit_s, bool show_rotation)
{
    // Initialise variables
    bool has_rotation = false;
    bool negative_rotation = false;
    // Loop over all positions, going row by row
    for (BDDVAR i = 0; i < circuit_s->nvars; i++) {
        for (BDDVAR j = 0; j < circuit_s->depth; j++) {
            // Print the gate
            printf("-%s",circuit_s->circuit[i][j].gateSymbol);
            // If rotation, check if a rotation needs to be printed in this column
            if (show_rotation) {
                // Reset rotation variables
                has_rotation = false;
                negative_rotation = false;
                // Check if any qubit in the column has a rotation and set variables accordingly
                for (BDDVAR k = 0; k < circuit_s->nvars; k++) {
                    if(circuit_s->circuit[k][j].id == gate_Rx.id || circuit_s->circuit[k][j].id == gate_Ry.id || circuit_s->circuit[k][j].id == gate_Rz.id)
                        has_rotation = true;
                    if(circuit_s->circuit[k][j].rotation < 0)
                        negative_rotation = true;
                }
                // If a rotation has been found in the column, print accordingly
                if (has_rotation) {
                    // If the current gate has a rotation, print its rotation
                    if(circuit_s->circuit[i][j].id == gate_Rx.id || circuit_s->circuit[i][j].id == gate_Ry.id || circuit_s->circuit[i][j].id == gate_Rz.id) {
                        printf("(%.4lf)",roundf(circuit_s->circuit[i][j].rotation*10000)/10000);
                        if(negative_rotation && circuit_s->circuit[i][j].rotation >= 0)
                            printf("-");
                    }
                    // If the current gate does not have a rotation, print extra wire space for alignment
                    else {
                        if (negative_rotation)
                            printf("---------");
                        else
                            printf("--------");
                    }
                }
            }
        }
        // After a wire has been printed, go to newline
        printf("-\n");
    }
}

void delete_circuit(Circuit* circuit_s)
{
    // Loop over all positions in the circuit
    for (BDDVAR j = 0; j < circuit_s->depth; j++) {
        for (BDDVAR i = 0; i < circuit_s->nvars; i++) {
            // If the gate has control qubits, free the list of indices
            if (circuit_s->circuit[i][j].id != gate_I.id) {
                if (circuit_s->circuit[i][j].control != NULL || circuit_s->circuit[i][j].controlSize != 0)
                    free(circuit_s->circuit[i][j].control);
            }
        }
    }
    // Free each wire
    for (BDDVAR i = 0; i < circuit_s->nvars; i++)
        free(circuit_s->circuit[i]);
    // Free the remaining data from the circuit struct
    free(circuit_s->circuit);
}

bool skip_to_gate(Circuit* circuit_s, BDDVAR i)
{
    // Initialise variables
    bool depth_reached, is_gate_I, is_gate_barrier;
    // Check if 'i' is valid
    if (i >= circuit_s->nvars)
        return false;
    // Walk over wire until the wire ends or you find a gate
    do {
        // Increase the progress of the wire
        circuit_s->progress[i]++;
        //Check contitions of next position and store results for while loop
        depth_reached = circuit_s->progress[i] < circuit_s->depth;
        is_gate_I = circuit_s->circuit[i][circuit_s->progress[i]].id == gate_I.id;
        is_gate_barrier = circuit_s->circuit[i][circuit_s->progress[i]].id == gate_barrier.id;
    } while (depth_reached && (is_gate_I || is_gate_barrier));
    return true;
}

TASK_IMPL_2(bool, advance, Circuit*, circuit_s, BDDVAR, q)
{
    // Initialise variables
    Gate gate = circuit_s->circuit[q][circuit_s->progress[q]];
    BDDVAR gate_id = get_gateid(gate);
    // Check if 'q' is valid
    if (q >= circuit_s->nvars)
        return false;
    // Check for every control qubit if the progress on that wire is equal to the progress of the target
    for (BDDVAR i = 0; i < gate.controlSize; i++) {
        if (circuit_s->circuit[gate.control[i]][circuit_s->progress[gate.control[i]]].id != gate_ctrl.id)
            return false;
    }
    // Handle single target and controlled gates
    if (gate.controlSize == 0)
        circuit_s->qdd = qdd_gate(circuit_s->qdd, gate_id, q);
    else if (gate.controlSize == 1)
        circuit_s->qdd = qdd_cgate(circuit_s->qdd, gate_id, gate.control[0], q);
    else if (gate.controlSize == 2)
        circuit_s->qdd = qdd_cgate2(circuit_s->qdd, gate_id, gate.control[0], gate.control[1], q);
    else if (gate.controlSize == 3)
        circuit_s->qdd = qdd_cgate3(circuit_s->qdd, gate_id, gate.control[0], gate.control[1], gate.control[2], q);
    // Advance to the successive gate on each wire that was affected above
    for (BDDVAR i = 0; i < gate.controlSize; i++)
        skip_to_gate(circuit_s, gate.control[i]);
    skip_to_gate(circuit_s, q);
    return true;
}

TASK_IMPL_1(BDDVAR, get_gateid, Gate, gate)
{
    // Initialise variables
    BDDVAR gate_id;
    // Switch case on the id of the gate_struct
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
        fprintf(stderr, "Unknown gate: %d\n", gate.id);
        exit(EXIT_FAILURE);
    }
    return gate_id;
}
