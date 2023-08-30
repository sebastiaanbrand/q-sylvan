#include <argp.h>
#include <sys/time.h>

#include "qsylvan.h"
#include "simple_parser.h"


/**********************<Arguments (configured via argp)>***********************/

static int workers = 1;
static int rseed = 0;
static bool count_nodes = false;
static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<16;
static size_t wgt_tab_size = 1LL<<23;
static double tolerance = 1e-14;
static int wgt_table_type = COMP_HASHMAP;
static int wgt_norm_strat = NORM_LARGEST;
static char* qasm_inputfile = NULL;
static char* json_outputfile = NULL;


static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed", 0},
    {"norm-strat", 's', "<low|largest|l2>", 0, "Edge weight normalization strategy", 0},
    {"tol", 't', "<tolerance>", 0, "Tolerance for deciding edge weights equal (default=1e-14)", 0},
    {"json", 'j', "<filename>", 0, "Write stats to given filename as json", 0},
    {"count-nodes", 'c', 0, 0, "Track maximum number of nodes", 1},
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
        else if (strcmp(arg, "largest")==0) wgt_norm_strat = NORM_LARGEST;
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


typedef struct stats_s {
    uint64_t applied_gates;
    uint64_t final_nodes;
    uint64_t max_nodes;
    uint64_t shots;
    double simulation_time;
} stats_t;
stats_t stats;


void fprint_stats(FILE *stream, quantum_circuit_t* circuit)
{
    fprintf(stream, "{\n");
    fprintf(stream, "  \"measurement_results\": {\n");
    fprintf(stream, "    \""); fprint_creg(stream, circuit); fprintf(stream, "\": 1\n");
    fprintf(stream, "  },\n");
    fprintf(stream, "  \"statistics\": {\n");
    fprintf(stream, "    \"applied_gates\": %ld,\n", stats.applied_gates);
    fprintf(stream, "    \"benchmark\": \"%s\",\n", circuit->name);
    fprintf(stream, "    \"final_nodes\": %ld,\n", stats.final_nodes);
    fprintf(stream, "    \"max_nodes\": %ld,\n", stats.max_nodes);
    fprintf(stream, "    \"n_qubits\": %d,\n", circuit->qreg_size);
    fprintf(stream, "    \"seed\": %d,\n", rseed);
    fprintf(stream, "    \"shots\": %ld,\n", stats.shots);
    fprintf(stream, "    \"simulation_time\": %lf,\n", stats.simulation_time);
    fprintf(stream, "    \"tolerance\": %.5e,\n", tolerance);
    fprintf(stream, "    \"wgt_norm_strat\": %d,\n", wgt_norm_strat);
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


QMDD apply_gate(QMDD state, quantum_op_t* gate)
{
    // This looks very ugly, but we need some way to match the name of a gate to
    // the GATEID, and since this is C we don't really have easy access to
    // data structures like a dictionary.
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
        return qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "cy") == 0) {
        return qmdd_cgate(state, GATEID_Y, gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        return qmdd_cgate(state, GATEID_Z, gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "ch") == 0) {
        return qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "csx") == 0) {
        return qmdd_cgate(state, GATEID_sqrtX, gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "crx") == 0) {
        return qmdd_cgate(state, GATEID_Rx(gate->angle[0]), gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "cry") == 0) {
        return qmdd_cgate(state, GATEID_Ry(gate->angle[0]), gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "crz") == 0) {
        return qmdd_cgate(state, GATEID_Rz(gate->angle[0]), gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "cp") == 0) {
        return qmdd_cgate(state, GATEID_Phase(gate->angle[0]), gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "cu") == 0) {
        return qmdd_cgate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->ctrls[0], gate->targets[0]);
    }
    else if (strcmp(gate->name, "ccx") == 0) {
        return qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
    }
    else if (strcmp(gate->name, "c3x") == 0) {
        return qmdd_cgate3(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0]);
    }
    else if (strcmp(gate->name, "c3sx") == 0) {
        return qmdd_cgate3(state, GATEID_sqrtX, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->targets[0]);
    }
    else if (strcmp(gate->name, "swap") == 0) {
        // no native SWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        return qmdd_circuit_swap(state, gate->targets[0], gate->targets[1]);
    }
    else if (strcmp(gate->name, "cswap") == 0) {
        // no native CSWAP gates in Q-Sylvan
        stats.applied_gates += 4;
        BDDVAR cs[3] = {gate->ctrls[0], AADD_INVALID_VAR, AADD_INVALID_VAR};
        return qmdd_ccircuit(state, CIRCID_swap, cs, gate->targets[0], gate->targets[1]);
    }
    else if (strcmp(gate->name, "rccx") == 0) {
        // no native RCCX (simplified Toffoli) gates in Q-Sylvan
        stats.applied_gates += 3;
        state = qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        state = qmdd_cgate2(state, GATEID_Z, gate->ctrls[0], gate->ctrls[1], gate->targets[0]);
        state = qmdd_gate(state, GATEID_X, gate->ctrls[1]);
        return state;
    }
    else if (strcmp(gate->name, "rzz") == 0 ) {
        // no native RZZ gates in Q-Sylvan
        stats.applied_gates += 2;
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        state = qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->targets[1]);
        state = qmdd_cgate(state, GATEID_X, gate->targets[0], gate->targets[1]);
        return state;
    }
    else if (strcmp(gate->name, "rxx") == 0) {
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
            state = apply_gate(state, op);
        }
        else if (op->type == op_measurement) {
            if (circuit->has_intermediate_measurements) {
                state = measure(state, op, circuit);
            }
            else {
                double p;
                state = qmdd_measure_all(state, circuit->qreg_size, circuit->creg, &p);
            }
        }
        if (count_nodes) {
            uint64_t count = aadd_countnodes(state);
            if (count > stats.max_nodes) stats.max_nodes = count;
        }
        op = op->next;
    }
    stats.simulation_time = wctime() - t_start;
    stats.shots = 1;
    stats.final_nodes = aadd_countnodes(state);
}


int main(int argc, char *argv[])
{
    argp_parse(&argp, argc, argv, 0, 0, 0);
    quantum_circuit_t* circuit = parse_qasm_file(qasm_inputfile);
    optimize_qubit_order(circuit);

    if (rseed == 0) rseed = time(NULL);
    srand(rseed);
    
    // Standard Lace initialization
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, tolerance, COMP_HASHMAP, wgt_norm_strat);

    simulate_circuit(circuit);

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
