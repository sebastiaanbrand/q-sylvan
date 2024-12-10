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

#include "qsylvan.h"
#include "qsylvan_qasm_parser.h"

/**********************<Arguments (configured via argp)>***********************/

static int workers = 1;
static int rseed = 0;
static bool count_nodes = false;
static bool output_vector = false;
static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<16;
static size_t min_wgt_tab_size = 1LL<<23;
static size_t max_wgt_tab_size = 1LL<<23;
static double tolerance = 1e-14;
static int wgt_table_type = COMP_HASHMAP;
static int wgt_norm_strat = NORM_MAX;
static bool wgt_inv_caching = true;
static int reorder_qubits = 0;
static char* qasm_inputfile = NULL;
static char* json_outputfile = NULL;


static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed", 0},
    {"norm-strat", 's', "<low|max|min|l2>", 0, "Edge weight normalization strategy", 0},
    {"tol", 't', "<tolerance>", 0, "Tolerance for deciding edge weights equal (default=1e-14)", 0},
    {"json", 'j', "<filename>", 0, "Write stats to given filename as json", 0},
    {"count-nodes", 'c', 0, 0, "Track maximum number of nodes", 0},
    {"state-vector", 'v', 0, 0, "Also output the complete state vector", 0},
    {"node-tab-size", 1000, "<size>", 0, "log2 of max node table size (max 40)", 0},
    {"wgt-tab-size", 1001, "<size>", 0, "log2 of max edge weigth table size (max 30 (23 if node table >2^30))", 0},
    {"reorder", 1002, 0, 0, "Reorders the qubits once such that (most) controls occur before targets in the variable order.", 0},
    {"reorder-swaps", 1003, 0, 0, "Reorders the qubits such that all controls occur before targets (requires inserting SWAP gates).", 0},
    {"disable-inv-caching", 1004, 0, 0, "Disable storing inverse of MUL and DIV in cache.", 0},
    {0, 0, 0, 0, 0, 0}
};

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
    case 's':
        if (strcmp(arg, "low")==0) wgt_norm_strat = NORM_LOW;
        else if (strcmp(arg, "max")==0) wgt_norm_strat = NORM_MAX;
        else if (strcmp(arg, "min")==0) wgt_norm_strat = NORM_MIN;
        else if (strcasecmp(arg, "l2")==0) wgt_norm_strat = NORM_L2;
        else argp_usage(state);
        break;
    case 't':
        tolerance = atof(arg);
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
    case 1000:
        if (atoi(arg) > 40) argp_usage(state);
        max_tablesize = 1LL<<(atoi(arg));
        break;
    case 1001:
        if (atoi(arg) > 30) argp_usage(state);
        max_wgt_tab_size = 1LL<<(atoi(arg));
        break;
    case 1002:
        reorder_qubits = 1;
        break;
    case 1003:
        reorder_qubits = 2;
        break;
    case 1004:
        wgt_inv_caching = false;
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

/*********************</Arguments (configured via argp)>***********************/


// NOTE: If there are only measurements at the end of the circuit, 'final_nodes' 
// and 'norm' will contain the node count and the norm of the state QMDD before
// the measurements.
typedef struct stats_s {
    uint64_t applied_gates;
    uint64_t final_nodes;
    uint64_t max_nodes;
    uint64_t shots;
    double simulation_time;
    double norm;
    QMDD final_state;
} stats_t;
stats_t stats;


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
            complex_t c = qmdd_get_amplitude(stats.final_state, x, circuit->qreg_size);
            fprintf(stream, "    [\n");
            fprintf(stream, "      %.16lf,\n", c.r);
            fprintf(stream, "      %.16lf\n", c.i);
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
    fprintf(stream, "    \"reorder\": %d,\n", reorder_qubits);
    fprintf(stream, "    \"seed\": %d,\n", rseed);
    fprintf(stream, "    \"shots\": %" PRIu64 ",\n", stats.shots);
    fprintf(stream, "    \"simulation_time\": %lf,\n", stats.simulation_time);
    fprintf(stream, "    \"tolerance\": %.5e,\n", tolerance);
    fprintf(stream, "    \"wgt_inv_caching\": %d,\n", wgt_inv_caching);
    fprintf(stream, "    \"wgt_norm_strat\": %d,\n", wgt_norm_strat);
    fprintf(stream, "    \"min_node_tab_size\": %" PRId64 ",\n", min_tablesize);
    fprintf(stream, "    \"max_node_tab_size\": %" PRId64 ",\n", max_tablesize);
    fprintf(stream, "    \"min_wgt_tab_size\": %" PRId64 ",\n", min_wgt_tab_size);
    fprintf(stream, "    \"max_wgt_tab_size\": %" PRId64 ",\n", max_wgt_tab_size);
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
 * Here we match the name of a gate in QASM to
 * the GATEID 
 */
QMDD apply_gate(QMDD state, quantum_op_t* gate, BDDVAR nqubits)
{
    // TODO: move this relation between parsed quantum_op and internal gate
    // somewhere else?
    stats.applied_gates++;

  if (strcmp(gate->name, "id") == 0) {
        stats.applied_gates--;
        return state;
    }
    else if (strcmp(gate->name, "x") == 0) {
        return qmdd_gate(state, GATEID_X, gate->targets[0]);
    }
    else if (strcmp(gate->name, "y") == 0) {
        return qmdd_gate(state, GATEID_Y, gate->targets[0]);
    }
    else if (strcmp(gate->name, "z") == 0) {
        return qmdd_gate(state, GATEID_Z, gate->targets[0]);
    }
    else if (strcmp(gate->name, "h") == 0) {
        return qmdd_gate(state, GATEID_H, gate->targets[0]);
    }
    else if (strcmp(gate->name, "s") == 0) {
        return qmdd_gate(state, GATEID_S, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sdg") == 0) {
        return qmdd_gate(state, GATEID_Sdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "t") == 0) {
        return qmdd_gate(state, GATEID_T, gate->targets[0]);
    }
    else if (strcmp(gate->name, "tdg") == 0) {
        return qmdd_gate(state, GATEID_Tdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sx") == 0) {
        return qmdd_gate(state, GATEID_sqrtX, gate->targets[0]);
    }
    else if (strcmp(gate->name, "sxdg") == 0) {
        return qmdd_gate(state, GATEID_sqrtXdag, gate->targets[0]);
    }
    else if (strcmp(gate->name, "rx") == 0) {
        return qmdd_gate(state, GATEID_Rx(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "ry") == 0) {
        return qmdd_gate(state, GATEID_Ry(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "rz") == 0) {
        return qmdd_gate(state, GATEID_Rz(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "p") == 0) {
        return qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "u2") == 0) {
        fl_t pi_over_2 = flt_acos(0.0);
        return qmdd_gate(state, GATEID_U(pi_over_2, gate->angle[0], gate->angle[1]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "u") == 0) {
        return qmdd_gate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->targets[0]);
    }
    else if (strcmp(gate->name, "cx") == 0) {
        return qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cy") == 0) {
        return qmdd_cgate(state, GATEID_Y, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        return qmdd_cgate(state, GATEID_Z, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "ch") == 0) {
        return qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "csx") == 0) {
        return qmdd_cgate(state, GATEID_sqrtX, gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "crx") == 0) {
        return qmdd_cgate(state, GATEID_Rx(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cry") == 0) {
        return qmdd_cgate(state, GATEID_Ry(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "crz") == 0) {
        return qmdd_cgate(state, GATEID_Rz(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cp") == 0) {
        return qmdd_cgate(state, GATEID_Phase(gate->angle[0]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "cu") == 0) {
        return qmdd_cgate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->ctrls[0], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "ccx") == 0) {
        return qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "c3x") == 0) {
        return qmdd_cgate3(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "c3sx") == 0) {
        return qmdd_cgate3(state, GATEID_sqrtX, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0], nqubits);
    }
    else if (strcmp(gate->name, "swap") == 0) {
        // no native SWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        return qmdd_circuit_swap(state, gate->targets[0], gate->targets[1]);
    }
    else if (strcmp(gate->name, "cswap") == 0) {
        // no native CSWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);
        // upside down CCNOT (equivalent)
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0], nqubits);
        // CCNOT
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->targets[0], gate->targets[1], nqubits);

        return state;
    }
    else if (strcmp(gate->name, "rccx") == 0) {
        // no native RCCX (simplified Toffoli) gates in Q-Sylvan
        stats.applied_gates += 3;
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->ctrls[1], gate->targets[0], nqubits);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        return state;
    }
    else if (strcmp(gate->name, "rzz") == 0 ) {
        // no native RZZ gates in Q-Sylvan
        stats.applied_gates += 2;
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);

        return state;
    }
    else if (strcmp(gate->name, "rxx") == 0) {
        // no native RXX gates in Q-Sylvan
        fl_t pi = flt_acos(0.0) * 2;
        stats.applied_gates += 6;
        state = qmdd_gate(state, GATEID_U(pi/2.0, gate->angle[0], 0), gate->targets[0]);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_Phase(-(gate->angle[0])), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1], nqubits);
        state = qmdd_gate(state, GATEID_H, gate->targets[1]);
        state = qmdd_gate(state, GATEID_U(pi/2.0, -pi, pi-gate->angle[0]), gate->targets[0]);
        return state;
    }
    else {
        fprintf(stderr, "Gate '%s' currently unsupported\n", gate->name);
        return state;
    }
}


QMDD measure(QMDD state, quantum_op_t *meas, quantum_circuit_t* circuit)
{
    double p;
    int m;
    printf("measure qubit %d, store result in creg[%d]\n", meas->targets[0], meas->meas_dest);
    qmdd_measure_qubit(state, meas->targets[0], circuit->qreg_size, &m, &p);
    circuit->creg[meas->meas_dest] = m;
    return state;
}


void simulate_circuit(quantum_circuit_t* circuit)
{
    double t_start = wctime();
    QMDD state = qmdd_create_all_zero_state(circuit->qreg_size);
    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        if (op->type == op_gate) {
            state = apply_gate(state, op, circuit->qreg_size);

        }
        else if (op->type == op_measurement) {
            if (circuit->has_intermediate_measurements) {
                state = measure(state, op, circuit);
            }
            else {
                double p;
                // don't set state = post measurement state
                qmdd_measure_all(state, circuit->qreg_size, circuit->creg, &p);
                if (circuit->reversed_qubit_order) {
                    reverse_bit_array(circuit->creg, circuit->qreg_size);
                }
                break;
            }
        }
        if (count_nodes) {
            uint64_t count = evbdd_countnodes(state);
            if (count > stats.max_nodes) stats.max_nodes = count;
        }
        op = op->next;
    }
    stats.simulation_time = wctime() - t_start;
    stats.final_state = state;
    stats.shots = 1;
    stats.final_nodes = evbdd_countnodes(state);
    stats.norm = qmdd_get_norm(state, circuit->qreg_size);
}


int main(int argc, char *argv[])
{
    argp_parse(&argp, argc, argv, 0, 0, 0);
    quantum_circuit_t* circuit = parse_qasm_file(qasm_inputfile);
    if (reorder_qubits)
        optimize_qubit_order(circuit, reorder_qubits == 2);

    if (rseed == 0) rseed = time(NULL);
    srand(rseed);
    
    // Standard Lace initialization
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(min_wgt_tab_size, max_wgt_tab_size, tolerance, COMP_HASHMAP, wgt_norm_strat);
    wgt_set_inverse_chaching(wgt_inv_caching);

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
