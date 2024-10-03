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

#include <inttypes.h>
#include <argp.h>
#include <sys/time.h>

#include <qsylvan.h>
#include <qsylvan_simulator_mtbdd.h>
#include <qsylvan_gates_mtbdd_mpc.h>

#include <sylvan_mpc.h>
#include "simple_parser.h" // TODO: rename in qsylvan_qasm_parser.h

/**
 * 
 * Arguments variables of command line interface (configured via argp).
 * 
 */
static int workers = 1;         // Number of threads running on separate CPU core
static int rseed = 0;
static int precision = 128;
static int rounding = 0;

static bool count_nodes = false;
static bool output_vector = false;

static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<16;

static char* qasm_inputfile = NULL;
static char* json_outputfile = NULL;

/**
 * 
 * Command Line Interface argument help list.
 * 
 */
static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed as integer", 0},
    {"precision", 'p', "<number of bits>", 0, "Precision of mantissa multiprecision complex float in bits (default=128)", 0},
    {"rounding", 'o', "<...>", 0, "Rounding strategy", 0},
    {"json", 'j', "<filename>", 0, "Write statistics in json format to <filename>", 0},
    {"count-nodes", 'c', 0, 0, "Track maximum number of nodes", 0},
    {"state-vector", 'v', 0, 0, "Show the complete state vector after simulation", 0},
    {0, 0, 0, 0, 0, 0}
};

/**
 * 
 * Command Line Interface argument processing.
 * 
 */

static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    switch (key) {

    case 'w':
        workers = atoi(arg);
        break;

    case 'r':
        rseed = atoi(arg);
        break;

    case 'p': // Precision
        precision = atoi(arg); // ASCII to integer, TODO: validate value
        break;

    case 'o': // Rounding
        rounding = atoi(arg);  // ASCII to integer, TODO: validate value
        break;

    case 'j':
        json_outputfile = arg;
        break;

    case 'c':
        count_nodes = true;
        break;

    case 'v':
        output_vector = true;
        break;

    case ARGP_KEY_ARG:
        if (state->arg_num >= 1) argp_usage(state);
        qasm_inputfile = arg;
        break;

    case ARGP_KEY_END:
        if (state->arg_num < 1) argp_usage(state);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
static struct argp argp = { options, parse_opt, "<qasm_file>", 0, 0, 0, 0 };


/**
 * 
 * Statistics of the simulation.
 * 
 * Store statistical info in struct below.
 * 
 * NOTE: If there are only measurements at the end of the circuit, 'final_nodes' 
 * and 'norm' will contain the node count and the norm of the state MTBDD before
 * the measurements.
 * 
 */

typedef struct stats_s {

    uint64_t applied_gates;     // Number of gates used
    uint64_t final_nodes;       // node count before the measurements
    uint64_t max_nodes;         // maximum of nodes >= final number of nodes
    uint64_t shots;             // Number of measurements
    double simulation_time;     // Time of simulation
    double norm;                // Norm L2 of the final state (should be 1.00000), TODO: make mpc type of this ?
    MTBDD final_state;          // State vector MTBDD after simulation

} stats_t;

stats_t stats;

/**
 * Print the statistics after the simulation the given file.
 */
void fprint_stats(FILE *stream, quantum_circuit_t* circuit)
{
    fprintf(stream, "{\n");
    fprintf(stream, "  \"measurement_results\": {\n");
    fprintf(stream, "    \""); fprint_creg(stream, circuit); fprintf(stream, "\": 1\n");
    fprintf(stream, "  },\n");
    if (output_vector)
    {
        fprintf(stream, "  \"state_vector\": [\n");

//printf("register size = %d\n", circuit->qreg_size);

        for (int k = 0; k < (1<<(circuit->qreg_size)); k++) {
            bool *x = int_to_bitarray(k, circuit->qreg_size, !(circuit->reversed_qubit_order));
            MTBDD leaf = mtbdd_getvalue_of_path(stats.final_state, x); // Perhaps reverse qubit sequence on circuit->qreg_size, see gmdd_get_amplitude
            fprintf(stream, "    [\n");
            mpc_out_str(stream, MPC_BASE_OF_FLOAT, 6, (mpc_ptr)mtbdd_getvalue(leaf), MPC_ROUNDING);
            if (k == (1<<(circuit->qreg_size))-1)
                fprintf(stream, "    ]\n");
            else
                fprintf(stream, "    ],\n");
            free(x);
        }
        fprintf(stream, "  ],\n");
    }
    fprintf(stream, "  \"statistics\": {\n");
    fprintf(stream, "    \"applied_gates\": %" PRIu64 ",\n", stats.applied_gates);
    fprintf(stream, "    \"benchmark\": \"%s\",\n", circuit->name);
    fprintf(stream, "    \"final_nodes\": %" PRIu64 ",\n", stats.final_nodes);
    fprintf(stream, "    \"max_nodes\": %" PRIu64 ",\n", stats.max_nodes);
    fprintf(stream, "    \"n_qubits\": %d,\n", circuit->qreg_size);
    fprintf(stream, "    \"norm\": %.5e,\n", stats.norm);
    fprintf(stream, "    \"seed\": %d,\n", rseed);
    fprintf(stream, "    \"shots\": %" PRIu64 ",\n", stats.shots);
    fprintf(stream, "    \"simulation_time\": %lf,\n", stats.simulation_time);
    fprintf(stream, "    \"precision\": %d,\n", precision);
    fprintf(stream, "    \"rounding\": %d,\n", rounding);
    fprintf(stream, "    \"workers\": %d\n", workers);
    fprintf(stream, "  }\n");
    fprintf(stream, "}\n");
}

static double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

/**
 * Here we match the QASM name of a gate with the MTBDD of the corresponding gate.
 * 
 * If the gate is supported the state vector |fi> will be applied to the matrix of the gate:
 * 
 *      M = (I (x) ... (x) I (x) G (x) I (x) ... (x) I) |fi> 
 * 
 * with I the 2x2 identity matrix, G the 2x2 choosen gate and |fi> the state vector.
 * 
 * The function returns the corresponding MTBDD of M.
 * 
 * The type of the matrix elements are complex numbers based on the GNU multiprecision complex library mpc.
 * 
 */
MTBDD apply_gate(MTBDD state, quantum_op_t* gate, int n)
{
    stats.applied_gates++;

    if (strcmp(gate->name, "id") == 0) {
        stats.applied_gates--;
        return state;
    }

    // Unitary none-parametrized gates

    else if (strcmp(gate->name, "x") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, X_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "y") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, Y_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "z") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, Z_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "h") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, H_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "s") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, S_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "sdg") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, S_dag_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "t") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, T_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "tdg") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, T_dag_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "sx") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, sqrt_X_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "sxdg") == 0) {
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, sqrt_X_dag_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    // Unitary parametrized gates

    else if (strcmp(gate->name, "rx") == 0) {
        MTBDD Rx_dd = mtbdd_Rx(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, Rx_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    
    else if (strcmp(gate->name, "ry") == 0) {
        MTBDD Ry_dd = mtbdd_Ry(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, Ry_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "rz") == 0) {
        MTBDD Rz_dd = mtbdd_Rz(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, Rz_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "p") == 0) {
        MTBDD P_dd = mtbdd_Phase(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    
//    else if (strcmp(gate->name, "u2") == 0) { // u3(a,b,c) = u(a,b,c) identical u2(b,c) = u(pi/2,b,c)
//        fl_t pi_over_2 = flt_acos(0.0);
//        return qmdd_gate(state, GATEID_U(pi_over_2, gate->angle[0], gate->angle[1]), gate->targets[0]);
//    }

    else if (strcmp(gate->name, "u") == 0) {
        MTBDD U_dd = mtbdd_U(gate->angle[0], gate->angle[1], gate->angle[2]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    // Control unitary none-parametrized gate combinations

    else if (strcmp(gate->name, "cx") == 0) {
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "cy") == 0) {
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, Y_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, Z_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "ch") == 0) {
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, H_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "csx") == 0) {
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    // Control unitary parametrized gate combinations

    else if (strcmp(gate->name, "crx") == 0) {
        MTBDD Rx_dd = mtbdd_Rx(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, Rx_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }    
    else if (strcmp(gate->name, "cry") == 0) {
        MTBDD Ry_dd = mtbdd_Ry(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, Ry_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "crz") == 0) {
        MTBDD Rz_dd = mtbdd_Rz(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, Rz_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "cp") == 0) {
        MTBDD P_dd = mtbdd_Phase(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, P_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "cu") == 0) {
        MTBDD U_dd = mtbdd_U(gate->angle[0], gate->angle[1], gate->angle[2]);
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, U_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }

    // Not used yet in benchmark circuits

/*
    else if (strcmp(gate->name, "ccx") == 0) {
        return qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
    }
    else if (strcmp(gate->name, "c3x") == 0) {
        return qmdd_cgate3(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0]);
    }
    else if (strcmp(gate->name, "c3sx") == 0) {
        return qmdd_cgate3(state, GATEID_sqrtX, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0]);
    }

    else if (strcmp(gate->name, "swap") == 0) { // swap(a,b) = cx(a,b); cx(b,a); cx(a,b)
        // no native SWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        return qmdd_circuit_swap(state, gate->targets[0], gate->targets[1]);
    }

    else if (strcmp(gate->name, "cswap") == 0) { // cswap(c,a,b) = cx(c,a,b); cx(c,b,a); cx(c,a,b)
        // no native CSWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1]);
        // upside down CCNOT (equivalent)
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0]);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->targets[0], gate->targets[1]);
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0]);
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1]);
        return state;
    }

    else if (strcmp(gate->name, "rccx") == 0) { // do not implement yet
        // no native RCCX (simplified Toffoli) gates in Q-Sylvan
        stats.applied_gates += 3;
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        return state;
    }
    else if (strcmp(gate->name, "rzz") == 0 ) { // do not implement yet
        // no native RZZ gates in Q-Sylvan
        stats.applied_gates += 2;
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        state = qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        return state;
    }
    else if (strcmp(gate->name, "rxx") == 0) { // do not implement yet
        // no native RXX gates in Q-Sylvan
        fl_t pi = flt_acos(0.0) * 2;
        stats.applied_gates += 6;
        state = qmdd_gate(state, GATEID_U(pi/2.0, gate->angle[0], 0), gate->targets[0]);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        state = qmdd_gate(state, GATEID_Phase(-(gate->angle[0])), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_gate(state, GATEID_U(pi/2.0, -pi, pi-gate->angle[0]), gate->targets[0]);
        return state;
    }
*/

    else {
        fprintf(stderr, "Gate '%s' currently unsupported\n", gate->name);
        return state;
    }
}

/**
 * 
 * Calculate the final state after observation (measurement).
 * 
 */
/* TODO: convert to MTBDD
MTBDD measure(QMDD state, quantum_op_t *meas, quantum_circuit_t* circuit)
{
    double p;
    int m;
    printf("measure qubit %d, store result in creg[%d]\n", meas->targets[0], meas->meas_dest);
    
// qmdd_measure_qubit(state, meas->targets[0], circuit->qreg_size, &m, &p);
    
    circuit->creg[meas->meas_dest] = m;
    return state;
}
*/

/**
 * 
 * Simulate the circuit as parsed.
 * 
 */
void simulate_circuit(quantum_circuit_t* circuit)
{
    double t_start = wctime();

    MTBDD state = mtbdd_create_all_zero_state_mpc(circuit->qreg_size);

    quantum_op_t *op = circuit->operations;
    while (op != NULL) {

        if (op->type == op_gate) {

            state = apply_gate(state, op, circuit->qreg_size); // , n);
        }

/* We only use the state vector and norm for now    
    
        else if (op->type == op_measurement) {

            assert(false && "Measurements not yet supported");
            fprintf(stderr, "Measurements not yet supported\n");
            exit(1);

            if (circuit->has_intermediate_measurements) {
                state = measure(state, op, circuit);
            }
            else {
                double p;
                // don't set state = post measurement state
                qmdd_measure_all(state, circuit->qreg_size, circuit->creg, &p); // TODO: convert to mtbdd 
                if (circuit->reversed_qubit_order) {
                    reverse_bit_array(circuit->creg, circuit->qreg_size);
                }
                break;
            }

        }
*/        
        if (count_nodes) {

uint64_t count = 0; // TODO: make mtbdd_countnodes(state);, exists.

            if (count > stats.max_nodes) stats.max_nodes = count;
        }
        op = op->next;
    }
    stats.simulation_time = wctime() - t_start;
    stats.final_state = state;
    stats.shots = 1;

//stats.final_nodes = mtbdd_countnodes(state);

    stats.norm = mtbdd_getnorm_mpc(state, circuit->qreg_size);

}

/**
 * 
 * Command Line Interface, parse arguments, parse QASM file, start simulation.
 * 
 */
int main(int argc, char *argv[])
{
    argp_parse(&argp, argc, argv, 0, 0, 0);

    //quantum_circuit_t* circuit = parse_qasm_file(qasm_inputfile);
    quantum_circuit_t* circuit = parse_qasm_file("/home/qrichard/Q-Sylvan/q-sylvan/qasm/circuits/test_circuit.qasm");

    print_quantum_circuit(circuit);

    optimize_qubit_order(circuit, false);

    if (rseed == 0) rseed = time(NULL);
    srand(rseed);
    
    // Standard Lace initialization
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    sylvan_init_mtbdd();

    // Set multi precision complex float type
    uint32_t mpc_type = mpc_init();
    assert(mpc_type == MPC_TYPE);

    // Init the dd's of the gates
    mtbdd_gates_init_mpc();

    printf("Simulate\n");
    simulate_circuit(circuit);
    print_quantum_circuit(circuit);

    printf("Statistics\n");
    fprint_stats(stdout, circuit);

    if (json_outputfile != NULL) {
        FILE *fp = fopen(json_outputfile, "w");
        fprint_stats(fp, circuit);
        fclose(fp);
    }

    sylvan_quit();
    lace_stop();
    free_quantum_circuit(circuit);

    return 0;
}
