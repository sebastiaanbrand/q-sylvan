#include "circuit_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"
#include <ctype.h>

/**
 * Prints the final results of the measured qubits.
 * 
 * @param qdd a QDD representing the statevector to be measured on
 * @param measurements an array containing booleans, where the to be measured qubits are flagged
 * @param nvars the number of qubits
 * @param runs the number of times the flagged qubits in the statevector should be measured
 * @param show toggle to print the measurement results
 */
void final_measure(QDD qdd, bool* measurements, BDDVAR nvars, BDDVAR runs, bool show)
{
    // Initialise variables
    bool *ms = malloc(nvars * sizeof(bool));
    double *p = malloc(nvars * sizeof(double));
    int *hits = malloc(pow(2,nvars) * sizeof(int));
    for(int i = 0; i < pow(2,nvars); i++) { hits[i] = 0; }

    // Run the circuit <runs> times
    for(BDDVAR i = 0; i < runs; i++) { 
        qdd_measure_all(qdd, nvars, ms, p);
        hits[bitarray_to_int(ms, nvars, false)]++;
    }

    // Reformat circuit results based on what qubits were measured
    BDDVAR index, j, sum = 0;
    for (BDDVAR i = 0; i < nvars; i++) { if (measurements[i]) sum++; }
    int *probs = malloc(pow(2, sum) * sizeof(int));
    for (int k = 0; k < pow(2, sum); k++) probs[k] = 0;
    for (int k = 0; k < (1 << (nvars)); k++) {
        bool *x = int_to_bitarray(k, nvars, false);
        j = sum-1;
        index = 0;
        for (BDDVAR i = nvars; i > 0; i--)
            if (measurements[i-1]) { index += x[i-1] * pow(2, j--); }
        probs[index] += hits[k];
        free(x);
    }

    // Print reformatted circuit results if <show> is toggled
    if (show) {
        for(int k = 0; k < pow(2, sum); k++) {
            if (probs[k] != 0) {
                bool *x = int_to_bitarray(k, nvars, false);
                j = sum-1;
                for (BDDVAR i = nvars; i > 0 ; i--) {
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

/**
 * Returns the Sylvan GATEID that corresponds to <gate>
 * 
 * @param gate Gate_struct of which to get the Sylvan GATEID
 * 
 * @return Sylvan GATEID which corresponds to <gate>
 */
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

/**
 * Applies <gate> to <qdd> on the qubit with index <i>.
 * 
 * @param qdd the statevector on which to apply <gate>
 * @param gate the gate to be applied
 * @param i the index of the qubit in <qdd> on which to apply <gate>
 * 
 * @return the statevector qdd where <gate> has been applied on the qubit with index <i>
 */
TASK_IMPL_3(QDD, apply_gate, QDD, qdd, Gate, gate, BDDVAR, i)
{
    // Get the corresponding Sylvan GATEID of <gate>
    BDDVAR gate_id = get_gate_id(gate);

    // Apply <gate> using the corresponding controls (if any)
    if (gate.controlSize == 0)
        qdd = qdd_gate(qdd, gate_id, i);
    else if (gate.controlSize == 1)
        qdd = qdd_cgate(qdd, gate_id, gate.control[0], i);
    else if (gate.controlSize == 2)
        qdd = qdd_cgate2(qdd, gate_id, gate.control[0], gate.control[1], i);
    else if (gate.controlSize == 3)
        qdd = qdd_cgate3(qdd, gate_id, gate.control[0], gate.control[1], gate.control[2], i);
    else {
        // More than 3 controls are not fully implemented yet, easier to throw on all cases
        fprintf(stderr, "Control-gates with more than 3 controls not implemented.\n");
        exit(EXIT_FAILURE);
    }
    return qdd;
}

/**
 * 
 */
TASK_IMPL_3(QDD, handle_control, Gate, gate, BDDVAR, k, BDDVAR, n)
{
    BDDVAR gate_id = get_gate_id(gate);
    int *c_options = malloc(n * sizeof(int));
    for (BDDVAR i = 0; i < n; i++) c_options[i] = -1;
    for (BDDVAR i = 0; i < gate.controlSize; i++) c_options[gate.control[i]] = 1;
    c_options[k] = 2;
    QDD qdd = qdd_create_multi_cgate_rec(n, c_options, gate_id, 0);
    return qdd;
}

QDD run_circuit_matrix(C_struct c_s, BDDVAR runs, bool show)
{
    LACE_ME;
    Gate gate;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    QDD vec = qdd_create_all_zero_state(c_s.nvars);
    QDD qdd, qdd_column;

    BDDVAR gate_id;
    BDDVAR *gateids = malloc(c_s.nvars * sizeof(BDDVAR));
    qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
    qdd_column = qdd;
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            gate = c_s.circuit[i][j];
            if (gate.id == gate_barrier.id)
                gateids[i] = GATEID_I;
            else if (gate.id == gate_measure.id) {
                gateids[i] = GATEID_I;
                measurements[i] = true;
            }
            else if (gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id)
                gateids[i] = GATEID_I;
            else if (gate.controlSize != 0) {
                gateids[i] = GATEID_I;
                qdd_column = handle_control(gate, i, c_s.nvars);
                qdd = qdd_matmat_mult(qdd_column, qdd, c_s.nvars);
            }
            else {
                gate_id = get_gate_id(gate);
                gateids[i] = gate_id;
            }
        }
        qdd_column = qdd_create_single_qubit_gates(c_s.nvars, gateids);
        if (gate.id != gate_barrier.id)
            qdd = qdd_matmat_mult(qdd_column, qdd, c_s.nvars);
        if (qdd_countnodes(qdd) > 100) {
            vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
            qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
        }
    }
    vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
    printf("%ld / %d\n", qdd_countnodes(vec), (int)pow(2,c_s.nvars));
    final_measure(vec, measurements, c_s.nvars, runs, show);
    return vec;
}

QDD run_c_struct(C_struct c_s, BDDVAR runs, bool show)
{
    LACE_ME;

    Gate gate;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    QDD qdd = qdd_create_all_zero_state(c_s.nvars);

    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            gate = c_s.circuit[i][j];
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_ctrl_c.id || gate.id == gate_I.id)
                continue;
            if (gate.id == gate_measure.id) {
                measurements[i] = true;
                continue;
            }
            qdd = apply_gate(qdd, gate, i);
        }
    }
    printf("%ld / %d\n", qdd_countnodes(qdd), (int)pow(2,c_s.nvars));
    final_measure(qdd, measurements, c_s.nvars, runs, show);
    free(measurements);
    // qdd_printnodes(qdd);
    return qdd;
}

/**
 * Runs QASM circuit given and prints the results
 * Input:
 *   the path to the file containing the QASM circuit code
 * Possible flags to be passed to the function:
 *   -r: the number of runs to perform
 *   -s: the randomness seed to be used
 *   -m: the boundaray value of nodes in a tree before multiplying with the state vector
 * Using the -m flag activates gate-gate multiplication runs, gate-statevector runs are used otherwise
 * 
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
    int opt;

    // Read flags from cmd and set parameters
    while((opt = getopt(argc, argv, "s:r:m:")) != -1) {
        switch(opt) {
            case 'r':
                if (isdigit(optarg))
                    runs = atoi(optarg);
                else
                    fprintf(stderr, "invalid argument for flag %d\n", opt);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'm':
                matrix = abs(atoi(optarg));
                break;
            default:
                fprintf(stderr, "usage: %s file [-r runs][-s seed][-m matrix_node_limit]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    printf("filename: %s\n", filename);
    printf("-r: %d\n", runs);
    printf("-s: %d\n", seed);
    printf("-m: %d\n", matrix);

    // Set randomness seed
    srand(seed);
    // Check if a file is given, if not, return an error
    if(access(filename, F_OK) != 0)
    {
        printf("Invalid QASM file.\n");
        exit(EXIT_FAILURE);
    }

    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<16, -1, COMP_HASHMAP);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    // Create a C_Struct representing the QASM circuit in the given file
    C_struct c_s = make_c_struct(filename, false);
    if (matrix != 0)
        run_circuit_matrix(c_s, runs, false);
    else
        run_c_struct(c_s, runs, false);

    delete_c_struct(&c_s);
    sylvan_quit();
    lace_exit();

    return 0;
}