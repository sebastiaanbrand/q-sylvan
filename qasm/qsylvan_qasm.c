#include <inttypes.h>
#include <popt.h>
#include <sys/time.h>
#include <libgen.h>

#include "qsylvan_qasm.h"

/**
 * Obtain current wallclock time
 */
static double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

static double t_start;
#define INFO(s, ...) fprintf(stdout, "[%8.2f] " s, wctime()-t_start, ##__VA_ARGS__)
#define Abort(...) { fprintf(stderr, __VA_ARGS__); fprintf(stderr, "Abort at line %d!\n", __LINE__); exit(-1); }

void final_measure(QMDD qmdd, int* measurements, C_struct c_s, bool* results)
{
    // Initialise variables
    bool *ms = malloc(c_s.qubits * sizeof(bool));
    double *p = malloc(c_s.qubits * sizeof(double));

    // Run the circuit <runs> times
    qmdd_measure_all(qmdd, c_s.qubits, ms, p);
    // Reformat circuit results based on what qubits were measured
    for (BDDVAR i = 0; i < c_s.qubits; i++) {
        if (measurements[i] != -1)
            results[measurements[i]] = ms[i];
    }

    // Free variables
    free(ms);
    free(p);
}

uint32_t get_gate_id(Gate gate)
{
    // Initialise variables
    BDDVAR gate_id;

    // Switch case on the id of the gate_struct
    if (gate.id == gate_I.id)
        gate_id = GATEID_I;
    else if (gate.id == gate_X.id)
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
    // Rx Ry and Rz gates have arbitrary rotations, Sylvan GATEIDs are allocated dynamically
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

QMDD apply_controlled_gate(QMDD qmdd, Gate gate, BDDVAR i, BDDVAR n)
{
    // Create controlled gate qmdd
    QMDD qmdd_column = handle_control_matrix(gate, i, n);
    // multiply with vector
    qmdd = aadd_matvec_mult(qmdd_column, qmdd, n);
    return qmdd;
}

QMDD apply_gate(QMDD qmdd, Gate gate, BDDVAR i, BDDVAR n)
{
    // Get the corresponding Sylvan GATEID of <gate>
    BDDVAR gate_id = get_gate_id(gate);
    // Apply <gate> using the corresponding controls (if any)
    if (gate.controlSize == 0)
        qmdd = qmdd_gate(qmdd, gate_id, i);
    // else if (gate.controlSize > 0 && gate.controlSize <= 4) {
    else
        qmdd = apply_controlled_gate(qmdd, gate, i, n);
    // }
    // else {
    //     // More than 3 controls are not fully implemented yet, easier to throw on all cases
    //     fprintf(stderr, "Control-gates with more than 3 controls not implemented.\n");
    //     exit(EXIT_FAILURE);
    // }
    return qmdd;
}

QMDD handle_control_matrix(Gate gate, BDDVAR k, BDDVAR n)
{
    // Initialise variables
    QMDD qmdd, qmdd_next;
    bool sorted = true;
    BDDVAR temp, gate_id, control = k, index = 0;
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    // Check if target is below controls
    for (BDDVAR j = 0; j < gate.controlSize; j++) {
        if (gate.control[j] > k)
            sorted = false;
        if (gate.control[j] > control) {
            control = gate.control[j];
            index = j;
        }
    }
    qmdd = qmdd_create_all_identity_matrix(n);
    // If target is not below some controls, swap target with lowest control
    if (!sorted) {
        qmdd = circuit_swap_matrix(k, control, n);
        temp = k;
        k = gate.control[index];
        gate.control[index] = temp;
    }
    // Get the corresponding Sylvan GATEID from <gate>
    gate_id = get_gate_id(gate);
    // Set the indexes of all controls of <gate> to 1
    for (BDDVAR i = 0; i < gate.controlSize; i++)
        c_options[gate.control[i]] = 1;
    // Set the index the target of <gate> to 2
    c_options[k] = 2;
    // Create the gate QMDD
    qmdd_next = qmdd_create_multi_cgate_rec(n, c_options, gate_id, 0);
    qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
    if (!sorted) {
        qmdd_next = circuit_swap_matrix(temp, k, n);
        qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
        temp = k;
        k = gate.control[index];
        gate.control[index] = temp;
    }
    return qmdd;
}

QMDD circuit_swap_matrix(BDDVAR qubit1, BDDVAR qubit2, BDDVAR n)
{
    // Initialise variables
    QMDD qmdd, qmdd_next;
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    // Set qubit1 to 1
    c_options[qubit1] = 1;
    // Set qubit2 to 2
    c_options[qubit2] = 2;

    // Perform swap
    qmdd = qmdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qmdd_next = qmdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
    qmdd_next = qmdd_create_multi_cgate_rec(n, c_options, GATEID_Z, 0);
    qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
    qmdd_next = qmdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
    qmdd_next = qmdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qmdd = aadd_matmat_mult(qmdd_next, qmdd, n);
    return qmdd;
}

bool is_final_measure(C_struct c_s, BDDVAR qubit, uint32_t depth)
{
    BDDVAR gateid;
    for (BDDVAR i = depth+1; i < c_s.depth; i++) {
        gateid = c_s.circuit[qubit][i].id;
        if (gateid != gate_barrier.id && gateid != gate_I.id) {
            return false;
        }
    }
    return true;
}

bool check_classical_if(uint32_t bits, Gate gate, bool* actual_list)
{
    int actual = 0;
    if (gate.classical_control == -1) {
        for (BDDVAR i = 0; i < bits; i++)
            actual += actual_list[i]*pow(2, i);
    }
    else
        actual = actual_list[gate.classical_control];
    return (gate.classical_expect == actual);
}

bool check_measuring_gates(C_struct c_s)
{
    for (BDDVAR i = 0; i < c_s.qubits; i++) {
        for (BDDVAR j = 0; j < c_s.depth; j++) {
            if (c_s.circuit[i][j].id == gate_measure.id) {
                for (BDDVAR k = j+1; k < c_s.depth; k++) {
                    if (c_s.circuit[i][k].id == gate_measure.id)
                        return true;
                }
            }
        }
    }
    return false;
}

void skip(C_struct c_s, BDDVAR* progress, BDDVAR i, bool skip_barrier)
{
    // Initialise variables
    bool depth_reached, is_gate_I, is_gate_ctrl_c, is_gate_barrier;
    // Walk over wire until the wire ends or you find a gate
    if (progress[i] > c_s.depth) return;
    do {
        // Increase the progress of the wire
        progress[i]++;
        //Check contitions of next position and store results for while loop
        depth_reached = progress[i] < c_s.depth;
        is_gate_I = c_s.circuit[i][progress[i]].id == gate_I.id;
        is_gate_ctrl_c = c_s.circuit[i][progress[i]].id == gate_ctrl_c.id;
        is_gate_barrier = c_s.circuit[i][progress[i]].id == gate_barrier.id;
    } while (depth_reached && (is_gate_I || is_gate_ctrl_c || (is_gate_barrier && skip_barrier)));
}

QMDD measure(QMDD qmdd, BDDVAR i, BDDVAR n, bool* result)
{
    int temp;
    double *p = malloc(sizeof(double));
    QMDD new_qmdd = qmdd_measure_qubit(qmdd, i, n, &temp, p);
    *result = (bool) temp;
    return new_qmdd;
}

bool check_gate(C_struct c_s, Gate gate, BDDVAR* progress, BDDVAR i, bool* results)
{
    bool ctrl_progress = true;
    if (gate.id == gate_ctrl.id || progress[i] >= c_s.depth)
        return false;
    // If some of its controls are blocked, the gate cannot be applied
    for (BDDVAR j = 0; j < gate.controlSize; j++) {
        if (progress[gate.control[j]] != progress[i])
            ctrl_progress = false;
    }
    if (!ctrl_progress)
        return false;
    // If classical control is not equal to max qubits, the gate is classically controlled
    if (gate.classical_expect != -1) {
        bool satisfied = check_classical_if(c_s.bits, gate, results);
        if (!satisfied)
            return false;
    }
    return true;
}

QMDD greedy(C_struct c_s, QMDD prev_qmdd, BDDVAR* column, int* measurements, bool* results, bool experiments, BDDVAR* n_gates)
{
    // Inisialise variables
    Gate gate, best_gate = gate_I;
    bool final, loop = true;
    BDDVAR best_nodecount, curr_nodecount, k = 0;
    QMDD best_qmdd = prev_qmdd, curr_qmdd = prev_qmdd;
    // Create a progress array
    BDDVAR* progress = malloc(c_s.qubits * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s.qubits; i++) progress[i] = *column;
    if (experiments) {
        best_nodecount = aadd_countnodes(prev_qmdd);
        printf("nodecount vec: %d at %d\n", best_nodecount, *n_gates);
    }
    for (BDDVAR i = 0; i < c_s.qubits; i++) {
        gate = c_s.circuit[i][progress[i]];
        if (gate.id == gate_I.id)
            skip(c_s, progress, i, false);
    }
    // Loop over all places in the circuit
    while(loop) {
        best_nodecount = pow(2, c_s.qubits+1);
        best_gate = gate_barrier;
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            gate = c_s.circuit[i][progress[i]];
            // if (gate.id == gate_I.id)
            //     skip(c_s, progress, i, false);
            if (!check_gate(c_s, gate, progress, i, results))
                continue;
            // Skip barrier and recheck qubit
            if (gate.id != gate_barrier.id) {
                // Set measurement flag
                if (gate.id == gate_measure.id) {
                    final = is_final_measure(c_s, i, progress[i]);
                    if (final)
                        measurements[i] = gate.control[0];
                    else {
                        curr_qmdd = measure(prev_qmdd, i, c_s.qubits, &results[gate.control[0]]);
                        // results[gate.control[0]] = (bool) result;
                    }
                }
                else
                    curr_qmdd = apply_gate(prev_qmdd, gate, i, c_s.qubits);
                curr_nodecount = aadd_countnodes(curr_qmdd);
                if (curr_nodecount < best_nodecount) {
                    best_nodecount = curr_nodecount;
                    best_qmdd = curr_qmdd;
                    best_gate = gate;
                    k = i;
                }
            }
        }
        if (best_gate.id != gate_barrier.id) {
            // continue to the next gate(s)
            if (best_gate.id != gate_measure.id) {
                for (BDDVAR i = 0; i < best_gate.controlSize; i++)
                    skip(c_s, progress, best_gate.control[i], false);
            }
            skip(c_s, progress, k, false);
            prev_qmdd = best_qmdd;
            *n_gates = *n_gates + 1;
            if (experiments)
                printf("nodecount vec: %d at %d\n", best_nodecount, *n_gates);
        }
        // Check if all wires are fully expanded (all gates are applied) or all qubits reached a barrier
        loop = false;
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            if (progress[i] < c_s.depth && c_s.circuit[i][progress[i]].id != gate_barrier.id)
                loop = true;
        }
    }
    best_nodecount = aadd_countnodes(prev_qmdd);
    if (experiments)
        printf("nodecount vec: %d at %d\n", best_nodecount, *n_gates);
    *column = progress[0];
    return prev_qmdd;
}

QMDD matmat(C_struct c_s, QMDD vec, BDDVAR* column, int* measurements, bool* results, int limit, bool experiments, BDDVAR* n_gates)
{
    Gate gate;
    bool final, satisfied;
    BDDVAR nodecount, j;
    BDDVAR *gateids = malloc(c_s.qubits * sizeof(BDDVAR));
    QMDD qmdd, qmdd_column;
    qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
    qmdd_column = qmdd;

    aadd_protect(&qmdd);
    aadd_protect(&qmdd_column);
    if (experiments) {
        nodecount = aadd_countnodes(qmdd);
        printf("nodecount mat: %d at %d\n", nodecount, *n_gates);
        nodecount = aadd_countnodes(vec);
        printf("nodecount vec: %d at %d\n", nodecount, *n_gates);
    }
    // Loop over all places in the circuit
    for (j = *column+1; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            gate = c_s.circuit[i][j];
            // Preset gateid (overwritten if needed)
            gateids[i] = GATEID_I;
            // Skip barrier (does not affect runs)
            // Skip control gates (controls are used by target gate)
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            // If classical control is not equal to -1, the gate is classically controlled
            if (gate.classical_expect != -1) {
                satisfied = check_classical_if(c_s.bits, gate, results);
                if (!satisfied)
                    continue;
            }
            // Set measurement flag
            if (gate.id == gate_measure.id) {
                final = is_final_measure(c_s, i, j);
                // If it is not a final measure, measure now and update QMDD
                if (final)
                    measurements[i] = gate.control[0];
                else {
                    // Multiply everything currently already in the column
                    qmdd_column = qmdd_create_single_qubit_gates(c_s.qubits, gateids);
                    qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
                    // Multiply with the vector and measure
                    vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
                    vec = measure(vec, i, c_s.qubits, &results[gate.control[0]]);
                    // Reset qmdd and gateids
                    qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
                    for(BDDVAR k = 0; k < i; k++)
                        gateids[k] = GATEID_I;
                }
            }
            // Handle controlled gates
            else if (gate.controlSize != 0) {
                qmdd_column = handle_control_matrix(gate, i, c_s.qubits);
                // Separately multiply with gate QMDD (not together with column)
                qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
            }
            // Set single gate in gateids
            else
                gateids[i] = get_gate_id(gate);
            *n_gates = *n_gates+1;
        }
        // If barrier, check if it starts a palindrome
        if (gate.id == gate_barrier.id)
            break;
        else {
            // Create gate QMDD out of current column (gateids)
            qmdd_column = qmdd_create_single_qubit_gates(c_s.qubits, gateids);
            // Multiply with previous columns
            qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
        }
        // Keep track of the all time highest node count in the matrices
        nodecount = aadd_countnodes(qmdd);
        // If nodecount of column QMDD is over limit...
        if (limit > 0 && nodecount > (BDDVAR)limit) {
            // Multiply with statevector QMDD and reset column QMDD
            vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
            qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
        }
        else if (experiments)
            printf("nodecount mat: %d at %d\n", nodecount, *n_gates);
    }
    *column = j+1;
    // Final multiply with statevector QMDD
    vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
    qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
    if (experiments) {
        nodecount = aadd_countnodes(qmdd);
        printf("nodecount mat: %d at %d\n", nodecount, *n_gates);
        nodecount = aadd_countnodes(vec);
        printf("nodecount vec: %d at %d\n", nodecount, *n_gates);
    }
    // Free variables
    free(gateids);
    aadd_unprotect(&qmdd);
    aadd_unprotect(&qmdd_column);
    return vec;
}

QMDD run_circuit_balance(C_struct c_s, int* measurements, bool* results, int limit, bool experiments)
{
    // Inisialise variables
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = -1;
    BDDVAR n_gates = 0, column = 0;
    QMDD vec = qmdd_create_all_zero_state(c_s.qubits);

    aadd_protect(&vec);
    while (column < c_s.depth) {
        vec = greedy(c_s, vec, &column, measurements, results, experiments, &n_gates);
        if (experiments)
            printf("palindrome start at: %d\n", n_gates);
        
        vec = matmat(c_s, vec, &column, measurements, results, limit, experiments, &n_gates);
        if (experiments && column < c_s.depth)
            printf("palindrome end at: %d\n", n_gates);
    }
    aadd_unprotect(&vec);
    return vec;
}

QMDD greedy_run_circuit(C_struct c_s, int* measurements, bool* results, bool experiments)
{
    // Initialize variables
    Gate gate, best_gate = gate_I;
    bool final, loop = true;
    BDDVAR best_nodecount, curr_nodecount, k = 0;
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = -1;
    QMDD prev_qmdd = qmdd_create_all_zero_state(c_s.qubits);
    QMDD best_qmdd = prev_qmdd, curr_qmdd = prev_qmdd;
    // Create a progress array
    BDDVAR* progress = malloc(c_s.qubits * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s.qubits; i++) progress[i] = 0;
    if (experiments) {
        best_nodecount = aadd_countnodes(prev_qmdd);
        printf("nodecount: %d\n", best_nodecount);
    }
    for (BDDVAR i = 0; i < c_s.qubits; i++) {
        gate = c_s.circuit[i][progress[i]];
        if (gate.id == gate_I.id)
            skip(c_s, progress, i, true);
    }
    // Loop over all places in the circuit
    while(loop) {
        best_nodecount = pow(2, c_s.qubits+1);
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            gate = c_s.circuit[i][progress[i]];
            if (!check_gate(c_s, gate, progress, i, results))
                continue;
            // Skip barrier and recheck qubit
            if (gate.id == gate_barrier.id) {
                skip(c_s, progress, i, true);
            }
            else {
                // Set measurement flag
                if (gate.id == gate_measure.id) {
                    final = is_final_measure(c_s, i, progress[i]);
                    if (final)
                        measurements[i] = gate.control[0];
                    else
                        curr_qmdd = measure(prev_qmdd, i, c_s.qubits, &results[gate.control[0]]);
                }
                else
                    curr_qmdd = apply_gate(prev_qmdd, gate, i, c_s.qubits);
                curr_nodecount = aadd_countnodes(curr_qmdd);
                if (curr_nodecount < best_nodecount) {
                    best_nodecount = curr_nodecount;
                    best_qmdd = curr_qmdd;
                    best_gate = gate;
                    k = i;
                }
            }
        }
        if (experiments)
            printf("nodecount: %d\n", best_nodecount);
        // continue to the next gate(s)
        if (best_gate.id != gate_measure.id) {
            for (BDDVAR i = 0; i < best_gate.controlSize; i++)
                skip(c_s, progress, best_gate.control[i], true);
        }
        skip(c_s, progress, k, true);
        prev_qmdd = best_qmdd;
        // Check if all wires are fully expanded (all gates are applied)
        loop = false;
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            if (progress[i] < c_s.depth)
                loop = true;
        }
    }
    return prev_qmdd;
}

QMDD run_circuit_matrix(C_struct c_s, int* measurements, bool* results, int limit, bool experiments)
{
    // Initialise variables
    // LACE_ME;
    Gate gate;
    bool final, satisfied, palindrome = false;
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = -1;
    BDDVAR gate_id, nodecount, n_gates = 0;
    BDDVAR *gateids = malloc(c_s.qubits * sizeof(BDDVAR));
    QMDD qmdd, qmdd_column;
    QMDD vec = qmdd_create_all_zero_state(c_s.qubits);
    qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
    qmdd_column = qmdd;

    aadd_protect(&vec);
    aadd_protect(&qmdd);
    aadd_protect(&qmdd_column);

    if (experiments) {
        nodecount = aadd_countnodes(qmdd);
        printf("nodecount mat: %d at %d\n", nodecount, n_gates);
        nodecount = aadd_countnodes(vec);
        printf("nodecount vec: %d at %d\n", nodecount, n_gates);
    }
    // Loop over all places in the circuit
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            gate = c_s.circuit[i][j];
            // Preset gateid (overwritten if needed)
            gateids[i] = GATEID_I;
            // Skip barrier (does not affect runs)
            // Skip control gates (controls are used by target gate)
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            // If classical control is not equal to -1, the gate is classically controlled
            if (gate.classical_expect != -1) {
                satisfied = check_classical_if(c_s.bits, gate, results);
                if (!satisfied)
                    continue;
            }
            // Set measurement flag
            if (gate.id == gate_measure.id) {
                final = is_final_measure(c_s, i, j);
                // If it is not a final measure, measure now and update QMDD
                if (final) {
                    measurements[i] = gate.control[0];
                }
                else {
                    // Multiply everything currently already in the column
                    qmdd_column = qmdd_create_single_qubit_gates(c_s.qubits, gateids);
                    qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
                    // Multiply with the vector and measure
                    vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
                    vec = measure(vec, i, c_s.qubits, &results[gate.control[0]]);
                    // Reset qmdd and gateids
                    qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
                    for(BDDVAR k = 0; k < i; k++)
                        gateids[k] = GATEID_I;
                }
            }
            // Handle controlled gates
            else if (gate.controlSize != 0) {
                qmdd_column = handle_control_matrix(gate, i, c_s.qubits);
                // Separately multiply with gate QMDD (not together with column)
                qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
            }
            // Set single gate in gateids
            else {
                gate_id = get_gate_id(gate);
                gateids[i] = gate_id;
            }
            n_gates++;
        }
        // If barrier, check if it starts a palindrome
        if (gate.id == gate_barrier.id) {
            palindrome = !palindrome;
            if (!palindrome && experiments)
                printf("palindrome end at: %d\n", n_gates);
            // get_palindrome_middle(c_s, j, &palindrome);
            if (palindrome && experiments)
                printf("palindrome start at: %d\n", n_gates);
        }
        else {
            // Create gate QMDD out of current column (gateids)
            qmdd_column = qmdd_create_single_qubit_gates(c_s.qubits, gateids);
            // Multiply with previous columns
            qmdd = aadd_matmat_mult(qmdd_column, qmdd, c_s.qubits);
        }
        // Keep track of the all time highest node count in the matrices
        // If nodecount of column QMDD is over limit...
        if (limit > 0) {
            nodecount = aadd_countnodes(qmdd);
            if (nodecount > (BDDVAR)limit) {
                // Multiply with statevector QMDD and reset column QMDD
                vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
                qmdd = qmdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
            }
        }
        if (experiments) {
            nodecount = aadd_countnodes(qmdd);
            printf("nodecount mat: %d at %d\n", nodecount, n_gates);
            nodecount = aadd_countnodes(vec);
            printf("nodecount vec: %d at %d\n", nodecount, n_gates);
        }
    }
    // Final multiply with statevector QMDD
    vec = aadd_matvec_mult(qmdd, vec, c_s.qubits);
    if (experiments) {
        nodecount = aadd_countnodes(vec);
        printf("nodecount vec: %d at %d\n", nodecount, n_gates);
    }
    // Free variables
    free(gateids);
    aadd_unprotect(&vec);
    aadd_unprotect(&qmdd);
    aadd_unprotect(&qmdd_column);
    return vec;
}

QMDD run_c_struct(C_struct c_s, int* measurements, bool* results, bool experiments)
{
    // Inisialise variables
    Gate gate;
    bool final, satisfied;
    double *p = malloc(sizeof(double));
    int result;
    BDDVAR nodecount;
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = -1;
    QMDD qmdd = qmdd_create_all_zero_state(c_s.qubits);
    aadd_protect(&qmdd);

    if (experiments) {
        nodecount = aadd_countnodes(qmdd);
        printf("nodecount: %d\n", nodecount);
    }
    // Loop over all places in the circuit
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.qubits; i++) {
            gate = c_s.circuit[i][j];
            // If classical control is not equal to max qubits, the gate is classically controlled
            if (gate.classical_expect != -1) {
                satisfied = check_classical_if(c_s.bits, gate, results);
                if (!satisfied)
                    continue;
            }
            // Skip barrier (does not affect runs)
            // Skip control gates (controls are used by target gate)
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            // Set measurement flag
            else if (gate.id == gate_measure.id) {
                final = is_final_measure(c_s, i, j);
                if (final)
                    measurements[i] = gate.control[0];
                else {
                    qmdd = qmdd_measure_qubit(qmdd, i, c_s.qubits, &result, p);
                    results[gate.control[0]] = (bool) result;
                }
            }
            // Apply gate
            else
                qmdd = apply_gate(qmdd, gate, i, c_s.qubits);
            if (experiments) {
                nodecount = aadd_countnodes(qmdd);
                printf("nodecount: %d\n", nodecount);
            }
        }
    }
    aadd_unprotect(&qmdd);
    return qmdd;
}

// TODO: move this main to separate file?
/**
 * Runs QASM circuit given by <filename> and prints the results.
 * Using the -m flag activates gate-gate multiplication runs, gate-statevector runs are used otherwise.
 * 
 * PARAMETERS:
 * - filename: the path to the file containing the QASM circuit code
 * 
 * FLAGS:
 * [-r runs (int)] (optional) the number of runs to perform
 * [-s seed (int)] (optional) the randomness seed to be used
 * [-m matrix (int)] (optional) the boundaray value of nodes in a tree before multiplying with the state vector
 * [-g greedy] (optional) runs the circuit matrix-vector method using a greedy algorithm
 * [-b balance (int)] (optional) runs the circuit switching between matrix-matrix method and greedy method
 * [-o optimize] (optional) optimize the circuit if true. This option will remove negating gates before running
 * [-e experiment] (optional) prints the nodcount and palindrome signals
 * [-t time] (optional) prints the time taken to run the circuit
 * 
 * NOTE:
 * Since multiplying a gate-QMDD with a gate-QMDD is more expensive than multiplying a gate-QMDD with
 * a statevector-QMDD, gate-gate multiplication is usually slower. However, circuits which use a lot of 
 * uncomputation can lead to smaller resulting gate QMDDs and possibly lead to faster runs.
 */
int main(int argc, char *argv[])
{
    t_start = wctime();
    // Initialise flag parameters
    int flag_help = -1;
    char *filename;
    int workers = 1;
    int wgt_norm_strat = NORM_LARGEST;
    unsigned int runs = 1;
    unsigned int seed = 0;
    int matrix = 0;
    int balance = 0;
    int greedy = 0;
    int optimize = 0;
    int intermediate_measuring = 0;
    int experiments = 0;
    bool intermediate_experiments;
    char *csv_outputfile = NULL;
    QMDD qmdd = AADD_TERMINAL;
    uint64_t res;

    poptContext con;
    struct poptOption optiontable[] = {
        { "help", 'h', POPT_ARG_NONE, &flag_help, 'h', "Display available options.", NULL },
        { "workers", 'w', POPT_ARG_INT, &workers, 'w', "Number of workers (cores). Default = 1.", NULL },
        { "runs", 'r', POPT_ARG_INT, &runs, 'r', "Number of runs to perform. Default = 1.", NULL },
        { "seed", 's', POPT_ARG_INT, &seed, 's', "Randomness seed to be used (!= 0). Default seeded with time().", NULL },
        { "matrix", 'm', POPT_ARG_INT, &matrix, 'm', "Boundaray value of nodes in a DD before multiplying with the state vector.", NULL },
        { "greedy", 'g', POPT_ARG_NONE, &greedy, 'g', "Runs the circuit matrix-vector method using a greedy algorithm.", NULL },
        { "balance", 'b', POPT_ARG_INT, &balance, 'b', "Runs the circuit switching between matrix-matrix method and greedy method", NULL },
        { "optimize", 'o', POPT_ARG_NONE, &optimize, 'o', "Optimize the circuit. This option will remove negating gates before running.", NULL },
        { "experiment", 'e', POPT_ARG_NONE, &experiments, 'e', "Prints the nodecount and palindrome signals.", NULL },
        { "norm-strat", 9, POPT_ARG_INT, &wgt_norm_strat, 9, "Weight norm strat as int: <0(low)|1(largest)|2(l2)>.", NULL },
        { "csv-output", 10, POPT_ARG_STRING, &csv_outputfile, 10, "Write stats to given filename (or append if file exists.", NULL },
        {NULL, 0, 0, NULL, 0, NULL, NULL}
    };
    con = poptGetContext("q-sylvan-sim", argc, (const char **)argv, optiontable, 0);
    poptSetOtherOptionHelp(con, "[OPTIONS..] <circuit.qasm>");
    
    if (argc < 2) {
        poptPrintUsage(con, stderr, 0);
        exit(1);
    }
    filename = argv[1];

    char c;  
    while ((c = poptGetNextOpt(con)) > 0) {
        switch (c) {
            case 'h':
                poptPrintHelp(con, stdout, 0);
                return 0;
        }
    }
    INFO("Option workers=%d\n", workers);
    INFO("Option norm-strat=%d\n", wgt_norm_strat);
    INFO("Option matrix=%d\n", matrix);
    INFO("Option balance=%d\n", balance);
    INFO("Option greedy=%d\n", greedy);
    INFO("Option experiments=%d\n", experiments);
    INFO("Option rseed=%d\n", seed);

    // Set randomness seed
    if (seed == 0)
        srand(time(NULL));
    else
         srand(seed);
    // Check if a file is given, if not, return an error
    if(access(filename, F_OK) != 0)
    {
        fprintf(stderr, "Invalid QASM file.\n");
        exit(EXIT_FAILURE);
    }

    double start,end;
    start = wctime();

    // Standard Lace initialization
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<23, 1LL<<23, -1, COMP_HASHMAP, wgt_norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests

    // Create a circuit struct representing the QASM circuit in the given file
    C_struct c_s = make_c_struct(filename, optimize);
    
    int* measurements = malloc(c_s.qubits * sizeof(int));
    bool* bit_res = malloc(c_s.bits * sizeof(bool));
    bool* bit_print;
    intermediate_measuring = check_measuring_gates(c_s);
    BDDVAR* results = malloc(pow(2,c_s.bits) * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < pow(2,c_s.bits); i++) results[i] = 0;
    BDDVAR* nodecount_matrix = malloc(sizeof(BDDVAR));
    *nodecount_matrix = 0;
    BDDVAR* nodecount = malloc(sizeof(BDDVAR));
    *nodecount = 0;

    if (intermediate_measuring) {
        intermediate_experiments = experiments;
        // Run the circuit based on method
        for (BDDVAR i = 0; i < runs; i++) {
            if (greedy)
                qmdd = greedy_run_circuit(c_s, measurements, bit_res, intermediate_experiments);
            else if (matrix != 0)
                qmdd = run_circuit_matrix(c_s, measurements, bit_res, matrix, intermediate_experiments);
            else if (balance != 0)
                qmdd = run_circuit_balance(c_s, measurements, bit_res, balance, intermediate_experiments);
            else
                qmdd = run_c_struct(c_s, measurements, bit_res, intermediate_experiments);
            // Measure all qubits
            final_measure(qmdd, measurements, c_s, bit_res);
            res = bitarray_to_int(bit_res, c_s.bits, false);
            results[res]++;
            experiments = false;
            intermediate_experiments = false;
        }
    }
    else {
        // Run the circuit based on method
        if (greedy)
                qmdd = greedy_run_circuit(c_s, measurements, bit_res, experiments);
        else if (matrix != 0)
            qmdd = run_circuit_matrix(c_s, measurements, bit_res, matrix, experiments);
        else if (balance != 0)
            qmdd = run_circuit_balance(c_s, measurements, bit_res, balance, experiments);
        else
            qmdd = run_c_struct(c_s, measurements, bit_res, experiments);
        // Measure all qubits
        for (BDDVAR i = 0; i < runs; i++) {
            final_measure(qmdd, measurements, c_s, bit_res);
            res = bitarray_to_int(bit_res, c_s.bits, false);
            results[res]++;
        }
    }

    // If experiments is true, print time
    end = wctime();
    double runtime = end-start;
    INFO("Circuit runtime: %lf\n", runtime);

    if (!experiments) {
        INFO("Measurement outcomes: (%d) measurements\n", runs);
        // Print reformatted circuit results if <show> is toggled
        for (BDDVAR i = 0; i < pow(2,c_s.bits); i++) {
            if (results[i] != 0) {
                bit_print = int_to_bitarray(i,c_s.bits,true);
                INFO("|");
                for (BDDVAR j = 0; j < c_s.bits; j++) { printf("%d", bit_print[j]); }
                printf("> : %d\n", results[i]);
            }
        }
    }

    uint64_t final_nodecount = aadd_countnodes(qmdd);
    double final_magnitude   = qmdd_get_magnitude(qmdd, c_s.qubits);
    INFO("Nodecount of final state: %" PRIu64 "\n", final_nodecount);
    INFO("Magnitude of final state: %.05lf\n", final_magnitude);
    
    if (csv_outputfile != NULL) {
        INFO("Writing stats to %s\n", csv_outputfile);

        FILE *fp = fopen(csv_outputfile, "a");
        // write header if file is empty
        fseek (fp, 0, SEEK_END);
            long size = ftell(fp);
            if (size == 0)
                fprintf(fp, "%s\n", "algorithm, nqubits, tolerance, norm-strat, inv-cache, workers, runtime, final_nodecount, final_magnitude");
        // append stats of this run
        char *alg_name = basename(filename);
        fprintf(fp, "%s, %d, %.3e, %d, %d, %d, %lf, %" PRIu64 ", %0.5lf\n",
                alg_name,
                c_s.qubits,
                wgt_store_get_tol(),
                weight_norm_strat,
                1, // inverse caching, TODO: make this configerable
                workers,
                runtime,
                final_nodecount,
                final_magnitude);
        fclose(fp);
        free(csv_outputfile);
    }

    // Free variables
    delete_c_struct(&c_s);
    sylvan_quit();
    lace_stop();

    return 0;
}
