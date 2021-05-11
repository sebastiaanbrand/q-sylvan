#include "circuit_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

void final_measure(QDD qdd, bool* measurements, BDDVAR qubits, BDDVAR runs, bool show)
{
    // Initialise variables
    bool *ms = malloc(qubits * sizeof(bool));
    double *p = malloc(qubits * sizeof(double));
    int *hits = malloc(pow(2,qubits) * sizeof(int));
    for(int i = 0; i < pow(2,qubits); i++) { hits[i] = 0; }

    // Run the circuit <runs> times
    for(BDDVAR i = 0; i < runs; i++) { 
        qdd_measure_all(qdd, qubits, ms, p);
        hits[bitarray_to_int(ms, qubits, false)]++;
    }

    // Reformat circuit results based on what qubits were measured
    BDDVAR index, j, sum = 0;
    for (BDDVAR i = 0; i < qubits; i++) { if (measurements[i]) sum++; }
    int *probs = malloc(pow(2, sum) * sizeof(int));
    for (int k = 0; k < pow(2, sum); k++) probs[k] = 0;
    for (int k = 0; k < (1 << (qubits)); k++) {
        bool *x = int_to_bitarray(k, qubits, false);
        j = sum-1;
        index = 0;
        for (BDDVAR i = qubits; i > 0; i--)
            if (measurements[i-1]) { index += x[i-1] * pow(2, j--); }
        probs[index] += hits[k];
        free(x);
    }

    // Print reformatted circuit results if <show> is toggled
    if (show) {
        for(int k = 0; k < pow(2, sum); k++) {
            if (probs[k] != 0) {
                bool *x = int_to_bitarray(k, qubits, false);
                j = sum-1;
                for (BDDVAR i = qubits; i > 0 ; i--) {
                    if (measurements[i-1]) {
                        printf("%d", x[j]);
                        j--;
                    }
                    else
                        printf("_");
                }
                printf(": %d\n", probs[k]);
                free(x);
            }
        }
    }

    // Free variables
    free(ms);
    free(p);
    free(hits);
    free(probs);
}

TASK_IMPL_1(BDDVAR, get_gate_id, Gate, gate)
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

TASK_IMPL_4(QDD, apply_controlled_gate, QDD, qdd, Gate, gate, BDDVAR, k, BDDVAR, n)
{
    // Create controlled gate qdd
    QDD qdd_column = handle_control_matrix(gate, k, n);
    // multiply with vector
    qdd = qdd_matvec_mult(qdd_column, qdd, n);
    return qdd;
}

TASK_IMPL_4(QDD, apply_gate, QDD, qdd, Gate, gate, BDDVAR, i, BDDVAR, n)
{
    // Get the corresponding Sylvan GATEID of <gate>
    BDDVAR gate_id = get_gate_id(gate);
    // Apply <gate> using the corresponding controls (if any)
    if (gate.controlSize == 0)
        qdd = qdd_gate(qdd, gate_id, i);
    // else if (gate.controlSize > 0 && gate.controlSize <= 4) {
    else
        qdd = apply_controlled_gate(qdd, gate, i, n);
    // }
    // else {
    //     // More than 3 controls are not fully implemented yet, easier to throw on all cases
    //     fprintf(stderr, "Control-gates with more than 3 controls not implemented.\n");
    //     exit(EXIT_FAILURE);
    // }
    return qdd;
}

TASK_IMPL_3(QDD, handle_control_matrix, Gate, gate, BDDVAR, k, BDDVAR, n)
{
    // Initialise variables
    QDD qdd, qdd_next;
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
    qdd = qdd_create_all_identity_matrix(n);
    // If target is not below some controls, swap target with lowest control
    if (!sorted) {
        qdd = circuit_swap_matrix(k, control, n);
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
    // Create the gate QDD
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, gate_id, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    if (!sorted) {
        qdd_next = circuit_swap_matrix(temp, k, n);
        qdd = qdd_matmat_mult(qdd_next, qdd, n);
        temp = k;
        k = gate.control[index];
        gate.control[index] = temp;
    }
    return qdd;
}

TASK_IMPL_3(QDD, circuit_swap_matrix, BDDVAR, qubit1, BDDVAR, qubit2, BDDVAR, n)
{
    // Initialise variables
    QDD qdd, qdd_next;
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    // Set qubit1 to 1
    c_options[qubit1] = 1;
    // Set qubit2 to 2
    c_options[qubit2] = 2;

    // Perform swap
    qdd = qdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qdd_next = qdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, GATEID_Z, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_single_qubit_gate(n, qubit1, GATEID_H);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    qdd_next = qdd_create_multi_cgate_rec(n, c_options, GATEID_X, 0);
    qdd = qdd_matmat_mult(qdd_next, qdd, n);
    return qdd;
}

TASK_IMPL_3(bool, is_final_measure, C_struct, c_s, BDDVAR, qubit, BDDVAR, depth)
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

TASK_IMPL_3(bool, check_classical_if, BDDVAR, bits, Gate, gate, int*, actual_list)
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

QDD run_circuit_matrix(C_struct c_s, BDDVAR runs, BDDVAR limit, bool show)
{
    // Inisialise variables
    LACE_ME;
    Gate gate;
    bool final, satisfied;
    double *p = malloc(sizeof(double));
    bool* measurements = malloc(c_s.qubits * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = false;
    int* results = malloc(c_s.bits * sizeof(int));
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    BDDVAR gate_id;
    BDDVAR *gateids = malloc(c_s.qubits * sizeof(BDDVAR));
    QDD qdd, qdd_column;
    QDD vec = qdd_create_all_zero_state(c_s.qubits);
    qdd = qdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
    qdd_column = qdd;

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
            // If classical control is not equal to max qubits, the gate is classically controlled
            if (gate.classical_expect != -1) {
                satisfied = check_classical_if(c_s.bits, gate, results);
                if (!satisfied)
                    continue;
            }
            // Set measurement flag
            if (gate.id == gate_measure.id) {
                final = is_final_measure(c_s, i, j);
                // If it is not a final measure, measure now and update QDD
                if (final) {
                    measurements[i] = true;
                }
                else {
                    // Multiply everything currently already in the column
                    qdd_column = qdd_create_single_qubit_gates(c_s.qubits, gateids);
                    qdd = qdd_matmat_mult(qdd_column, qdd, c_s.qubits);
                    // Multiply with the vector and measure
                    vec = qdd_matvec_mult(qdd, vec, c_s.qubits);
                    vec = qdd_measure_qubit(vec, i, c_s.qubits, &results[gate.control[0]], p);
                    // Reset qdd and gateids
                    qdd = qdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
                    for(BDDVAR k = 0; k < i; k++)
                        gateids[k] = GATEID_I;
                }
            }
            // Handle controlled gates
            else if (gate.controlSize != 0) {
                qdd_column = handle_control_matrix(gate, i, c_s.qubits);
                // Separately multiply with gate QDD (not together with column)
                qdd = qdd_matmat_mult(qdd_column, qdd, c_s.qubits);
            }
            // Set single gate in gateids
            else {
                gate_id = get_gate_id(gate);
                gateids[i] = gate_id;
            }
        }
        // Skip barrier (whole column has no effect)
        if (gate.id != gate_barrier.id) {
            // Create gate QDD out of current column (gateids)
            qdd_column = qdd_create_single_qubit_gates(c_s.qubits, gateids);
            // Multiply with previous columns
            qdd = qdd_matmat_mult(qdd_column, qdd, c_s.qubits);
        }
        // If nodecount of column QDD is over limit...
        if (qdd_countnodes(qdd) > limit) {
            // Multiply with statevector QDD and reset column QDD
            vec = qdd_matvec_mult(qdd, vec, c_s.qubits);
            qdd = qdd_create_single_qubit_gates_same(c_s.qubits, GATEID_I);
        }
    }
    // Final multiply with statevector QDD
    vec = qdd_matvec_mult(qdd, vec, c_s.qubits);
    // Measure all qubits
    final_measure(vec, measurements, c_s.qubits, runs, show);
    // Free variables
    free(measurements);
    free(gateids);
    return vec;
}

QDD run_c_struct(C_struct c_s, BDDVAR runs, bool show)
{
    // Inisialise variables
    LACE_ME;
    Gate gate;
    bool final, satisfied;
    double *p = malloc(sizeof(double));
    bool* measurements = malloc(c_s.qubits * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.qubits; i++) measurements[i] = false;
    int* results = malloc(c_s.bits * sizeof(int));
    for (BDDVAR i = 0; i < c_s.bits; i++) results[i] = 0;
    QDD qdd = qdd_create_all_zero_state(c_s.qubits);

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
                    measurements[i] = true;
                else
                    qdd = qdd_measure_qubit(qdd, i, c_s.qubits, &results[gate.control[0]], p);
            }
            // Apply gate
            else
                qdd = apply_gate(qdd, gate, i, c_s.qubits);
        }
    }
    // Measure all qubits
    final_measure(qdd, measurements, c_s.qubits, runs, show);
    // Free variables
    free(measurements);
    return qdd;
}

/**
 * Runs QASM circuit given by <filename> and prints the results.
 * Using the -m flag activates gate-gate multiplication runs, gate-statevector runs are used otherwise.
 * 
 * PARAMETERS:
 * - filename: the path to the file containing the QASM circuit code
 * 
 * FLAGS:
 * [-r runs] (optional) the number of runs to perform
 * [-s seed] (optional) the randomness seed to be used
 * [-m matrix] (optional) the boundaray value of nodes in a tree before multiplying with the state vector
 * [-o optimize] (optional) optimize the circuit if true. This option will remove negating gates before running
 * 
 * NOTE:
 * Since multiplying a gate-QDD with a gate-QDD is more expensive than multiplying a gate-QDD with
 * a statevector-QDD, gate-gate multiplication is usually slower. However, circuits which use a lot of 
 * uncomputation can lead to smaller resulting gate QDDs and possibly lead to faster runs.
 */
int main(int argc, char *argv[])
{
    // Initialise flag parameters
    char *filename = argv[1];
    BDDVAR runs = 100;
    BDDVAR seed = 100;
    BDDVAR matrix = 0;
    bool optimize = false;
    int opt;

    // Read flags from cmd and set parameters
    while((opt = getopt(argc, argv, "s:r:m:o")) != -1) {
        switch(opt) {
            case 'r':
                runs = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'm':
                matrix = abs(atoi(optarg));
                break;
            case 'o':
                optimize = true;
                break;
            default:
                fprintf(stderr, "usage: %s file [-r runs][-s seed][-m matrix_node_limit][-o optimize]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Set randomness seed
    srand(seed);
    // Check if a file is given, if not, return an error
    if(access(filename, F_OK) != 0)
    {
        fprintf(stderr, "Invalid QASM file.\n");
        exit(EXIT_FAILURE);
    }

    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<16, -1, COMP_HASHMAP, NORM_LARGEST);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    // Create a circuit struct representing the QASM circuit in the given file
    C_struct c_s = make_c_struct(filename, optimize);
    // Run the circuit based on method
    if (matrix != 0)
        run_circuit_matrix(c_s, runs, matrix, true);
    else
        run_c_struct(c_s, runs, true);
    // Free variables
    delete_c_struct(&c_s);
    sylvan_quit();
    lace_exit();

    return 0;
}