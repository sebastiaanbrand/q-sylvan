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
#include "qsylvan_qasm_parser.h"

/**
 * 
 * Arguments variables of command line interface (configured via argp).
 * 
 */
static int workers = 1;         // Number of threads running on separate CPU core
static int rseed = 0;
static int precision = MPC_PRECISION;
static int rounding = 0;

static bool count_nodes = true;
static bool output_vector = false;

static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<16;

static char* qasm_inputfile = NULL;
static char* json_outputfile = NULL;

#define M_PI 3.14159265358979323846
static double PI_1 = M_PI;
static double PI_2 = M_PI / 2.0;
static double PI_4 = M_PI / 4.0;
static double PI_8 = M_PI / 8.0;

/**
 * 
 * Command Line Interface argument help list.
 * 
 */
static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed as integer", 0},
    {"precision", 'p', "<number of bits>", 0, "Precision of mantissa multiprecision complex float in bits (default=MPC_PRECISION)", 0},
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

        for (int k = 0; k < (1<<(circuit->qreg_size)); k++) {
            bool *x = int_to_bitarray(k, circuit->qreg_size, !(circuit->reversed_qubit_order));
            MTBDD leaf = mtbdd_getvalue_of_path(stats.final_state, x); // Perhaps reverse qubit sequence on circuit->qreg_size, see gmdd_get_amplitude
            fprintf(stream, "    [\n");
            
            mpfr_t real, imag;
            mpfr_init2(real, MPC_PRECISION);
            mpfr_init2(imag, MPC_PRECISION);
            mpc_real(real, (mpc_ptr)mtbdd_getvalue(leaf), MPC_ROUNDING);
            mpc_imag(imag, (mpc_ptr)mtbdd_getvalue(leaf), MPC_ROUNDING);
            mpfr_fprintf(stream, "      %.17Rf,\n", real);
            mpfr_fprintf(stream, "      %.17Rf\n", imag);
            //mpc_out_str(stream, MPC_BASE_OF_FLOAT, 17, (mpc_ptr)mtbdd_getvalue(leaf), MPC_ROUNDING);
            
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
    else if (strcmp(gate->name, "u1") == 0) {
        MTBDD U_dd = mtbdd_U1(gate->angle[0]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "u2") == 0) {
        MTBDD U_dd = mtbdd_U2(gate->angle[0], gate->angle[1]);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        return mtbdd_matvec_mult(M_dd, state, 2*n, 0);
    }
    else if (strcmp(gate->name, "u3") == 0 || strcmp(gate->name, "u") == 0) {
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

    // Composed gates based on (control) unitary gates

    else if (strcmp(gate->name, "ccx") == 0) {
        MTBDD CCX_dd = mtbdd_ccg(n, gate->ctrls[0], gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        return mtbdd_matvec_mult(CCX_dd, state, 2*n, 0);
    }

    else if (strcmp(gate->name, "c3x") == 0) {

        // c3x a,b,c,d = c3x ctrls[0], ctrls[1], ctrls[2], targets[0]

        // h d;
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) a;
        MTBDD P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) b;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[1], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) c;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) b;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[1], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->ctrls[2], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) c;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[2], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) c;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->ctrls[2], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) c;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[2], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx b, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx b, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // h d;
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 31 -  1;

        return state;
    }
    else if (strcmp(gate->name, "c3sx") == 0) {

        // c3sx a,b,c,d = c3x ctrls[0], ctrls[1], ctrls[2], targets[0]

        // h d;
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) a;
        MTBDD P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) b;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[1], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) c;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[1], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) b;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[1], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[1], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->ctrls[2], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) c;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[2], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) c;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->ctrls[2], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) c;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->ctrls[2], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->ctrls[2], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx b, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx b, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(pi/8) d;
        P_dd = mtbdd_Phase(PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx c, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[2], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // p(-pi/8) d;
        P_dd = mtbdd_Phase(-PI_8);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, P_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // csx a, d;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, sqrt_X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // h d;
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 31 - 1;

        return state;
    }
    else if (strcmp(gate->name, "swap") == 0) {

        // swap(a,b) = cx(a,b); cx(b,a); cx(a,b), swap(|q0> (x) |q1>) = |q1> (x) |q0>

        MTBDD M1_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M1_dd, state, 2*n, 0);

        MTBDD M2_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M2_dd, state, 2*n, 0);

        MTBDD M3_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);

        state = mtbdd_matvec_mult(M3_dd, state, 2*n, 0); 

        stats.applied_gates += 3 - 1;

        return state;
    }
    else if (strcmp(gate->name, "cswap") == 0) { 

        // cswap a,b,c = csawp ctrls[0], target[0], target[1]

        // cx c,b;
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // ccx a,b,c;
        MTBDD CCX_dd = mtbdd_ccg(n, gate->ctrls[0], gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(CCX_dd, state, 2*n, 0);
        
        // cx c,b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 3 - 1;

        return state;
    }
    else if (strcmp(gate->name, "rccx") == 0) {

        // rccx a,b,c = rccx ctrls[0] ctrls[1] targets[0]

        // u2(0,pi) c;
        MTBDD U_dd = mtbdd_U2(0.0, PI_1);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
  
        // u1(pi/4) c;
        U_dd = mtbdd_U1(PI_4);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
  
        // cx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // u1(-pi/4) c;
        U_dd = mtbdd_U1(-PI_4);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
  
        // cx a, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[0], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
  
        // u1(pi/4) c;
        U_dd = mtbdd_U1(PI_4);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx b, c;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->ctrls[1], gate->targets[0], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // u1(-pi/4) c;
        U_dd = mtbdd_U1(-PI_4);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // u2(0,pi) c;
        U_dd = mtbdd_U2(0.0, PI_1);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 9 - 1;

        return state;
    }
    else if (strcmp(gate->name, "rzz") == 0 ) {
    
        // rzz(theta) a,b = rzz(theta) targets[0] targets[1]

        // cx a,b;
        MTBDD M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // u1(theta) b;
        MTBDD U_dd = mtbdd_U1(gate->angle[0]);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[1], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // cx a,b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 3 - 1;

        return state;
    }
    else if (strcmp(gate->name, "rxx") == 0) {

        // rxx(theta) a,b

        // u3(pi/2, theta, 0) a;
        MTBDD U_dd = mtbdd_U(PI_2, gate->angle[0], 0.0);
        MTBDD M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // h b;
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[1], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // cx a,b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // u1(-theta) b;
        U_dd = mtbdd_U1(-gate->angle[0]);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        // cx a,b;
        M_dd = mtbdd_create_single_control_gate_for_qubits_mpc(n, gate->targets[0], gate->targets[1], I_dd, V00_dd, V11_dd, X_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // h b;
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[1], I_dd, H_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);
        
        // u2(-pi, pi-theta) a;
        U_dd = mtbdd_U2(-PI_1, PI_1 - gate->angle[0]);
        M_dd = mtbdd_create_single_gate_for_qubits_mpc(n, gate->targets[0], I_dd, U_dd);
        state = mtbdd_matvec_mult(M_dd, state, 2*n, 0);

        stats.applied_gates += 7 - 1;

        return state;
    }

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
    
    qmdd_measure_qubit(state, meas->targets[0], circuit->qreg_size, &m, &p);
    
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

            uint64_t count = mtbdd_nodecount(state);
            if (count > stats.max_nodes) stats.max_nodes = count;
        }
        op = op->next;
    }
    stats.simulation_time = wctime() - t_start;
    stats.final_state = state;
    stats.shots = 1;

    stats.final_nodes = mtbdd_nodecount(state);
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

    quantum_circuit_t* circuit = parse_qasm_file(qasm_inputfile);

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

    simulate_circuit(circuit);

    if (json_outputfile != NULL) {
        FILE *fp = fopen(json_outputfile, "w");
        fprint_stats(fp, circuit);
        fclose(fp);
    } else {
        fprint_stats(stdout, circuit);
    }

    sylvan_quit();
    lace_stop();
    free_quantum_circuit(circuit);

    return 0;
}
