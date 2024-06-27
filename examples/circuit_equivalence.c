#include <inttypes.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <time.h>
#include <sys/time.h>
#include <qsylvan.h>
#include "../qasm/simple_parser.h"

#define max(a,b) ((a > b) ? a : b)

static int workers = 1;
static size_t min_tablesize = 1LL<<18;
static size_t max_tablesize = 1LL<<20;
static size_t min_cachesize = 1LL<<16;
static size_t max_cachesize = 1LL<<18;
static size_t wgt_tab_size = 1LL<<15;
static double tolerance = 1e-14;
static int wgt_table_type = COMP_HASHMAP;
static int wgt_norm_strat = NORM_MAX;
static int count_nodes = 1;


typedef struct stats_s {
    bool are_equiv;
    double cpu_time;
    double wall_time;
    size_t max_nodes_total;
    size_t max_nodes_U;
    size_t max_nodes_V;
} stats_t;
stats_t stats = {0};

quantum_circuit_t* circuit_U;
quantum_circuit_t* circuit_V;


void print_stats() {
    // print stats in JSON format
    printf("{\n");
    printf("  \"statistics\": {\n");
    printf("    \"are_equiv\" : %d,\n", (int)stats.are_equiv);
    printf("    \"circuit_U\": \"%s\",\n", circuit_U->name);
    printf("    \"circuit_V\": \"%s\",\n", circuit_V->name);
    printf("    \"max_nodes_total\": %" PRIu64 ",\n", stats.max_nodes_total);
    printf("    \"max_nodes_U\": %" PRIu64 ",\n", stats.max_nodes_U);
    printf("    \"max_nodes_V\": %" PRIu64 ",\n", stats.max_nodes_V);
    printf("    \"n_qubits\": %d,\n", circuit_U->qreg_size);
    printf("    \"time_cpu\": %lf,\n", stats.cpu_time);
    printf("    \"time_wall\": %lf,\n", stats.wall_time);
    printf("    \"tolerance\": %.5e,\n", tolerance);
    printf("    \"wgt_norm_strat\": %d,\n", wgt_norm_strat);
    printf("    \"workers\": %d\n", workers);
    printf("  }\n");
    printf("}\n");
}

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


QMDD get_gate_matrix(quantum_op_t* gate, BDDVAR nqubits, bool dag) {
    if (strcmp(gate->name, "id") == 0) {
        return qmdd_create_all_identity_matrix(nqubits);
    }
    else if (strcmp(gate->name, "x") == 0) {
        return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_X);
    }
    else if (strcmp(gate->name, "y") == 0) {
        return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Y);
    }
    else if (strcmp(gate->name, "z") == 0) {
        return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Z);
    }
    else if (strcmp(gate->name, "h") == 0) {
        return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_H);
    }
    else if (strcmp(gate->name, "s") == 0) {
        if (dag) return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Sdag);
        else     return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_S);
    }
    else if (strcmp(gate->name, "sdg") == 0) {
        if (dag) return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_S);
        else     return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Sdag);
    }
    else if (strcmp(gate->name, "t") == 0) {
        if (dag) return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Tdag);
        else     return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_T);
    }
    else if (strcmp(gate->name, "tdg") == 0) {
        if (dag) return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_T);
        else     return qmdd_create_single_qubit_gate(nqubits, gate->targets[0], GATEID_Tdag);
    }
    else if (strcmp(gate->name, "cx") == 0) {
        return qmdd_create_controlled_gate(nqubits, gate->ctrls[0], gate->targets[0], GATEID_X);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        return qmdd_create_controlled_gate(nqubits, gate->ctrls[0], gate->targets[0], GATEID_Z);
    }
    else {
        fprintf(stderr, "Gate '%s' currently unsupported\n", gate->name);
        exit(EXIT_FAILURE);
    }
}

QMDD compute_UPUdag(quantum_circuit_t *circuit, gate_id_t P, BDDVAR k) {
    BDDVAR nqubits = circuit->qreg_size;
    QMDD circ_matrix = qmdd_create_single_qubit_gate(nqubits, k, P);
    QMDD tmp = AADD_ZERO;
    aadd_protect(&circ_matrix);
    aadd_protect(&tmp);

    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        switch (op->type) {
            case op_gate:
                tmp = get_gate_matrix(op, nqubits, false);
                circ_matrix = aadd_matmat_mult(tmp, circ_matrix, nqubits);
                tmp = get_gate_matrix(op, nqubits, true);
                circ_matrix = aadd_matmat_mult(circ_matrix, tmp, nqubits);
                break;
            case op_measurement:
                fprintf(stderr, "ERROR: Measurments not supported\n");
                exit(EXIT_FAILURE);
                break;
            default:
                break;
        }
        op = op->next;
    }

    aadd_unprotect(&circ_matrix);
    aadd_unprotect(&tmp);

    return circ_matrix;
}

 void circuit_equivalence_check(quantum_circuit_t *U, quantum_circuit_t *V) {
    if (U->qreg_size != V->qreg_size) {
        fprintf(stderr, "Circuits have different number of qubits\n");
        return;
    }

    BDDVAR nqubits = U->qreg_size;

    uint64_t m = 0;
    uint32_t XZ[2] = {GATEID_X, GATEID_Z};
    char XZ_str[2] = {'X', 'Z'};
    for (BDDVAR k = 0; k < nqubits; k++) {
        for (int i = 0; i < 2; i++) {
            QMDD qmdd_U = compute_UPUdag(U, XZ[i], k);
            QMDD qmdd_V = compute_UPUdag(V, XZ[i], k);

            if (count_nodes) {
                size_t nodes_U = aadd_countnodes(qmdd_U);
                size_t nodes_V = aadd_countnodes(qmdd_V);
                stats.max_nodes_total = max(stats.max_nodes_total,nodes_U+nodes_V);
                stats.max_nodes_U = max(stats.max_nodes_U, nodes_U);
                stats.max_nodes_V = max(stats.max_nodes_V, nodes_V);
            }

            if (qmdd_U != qmdd_V) {
                stats.are_equiv = false;
                printf("U*%c_%d*U^dag != V*%c_%d*V^dag (%zu nodes)\n",
                        XZ_str[i], k, XZ_str[i], k, (size_t)m);
                return;
            }
        }
    }

    stats.are_equiv = true;
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Expected ./circuit_eq file1.qasm file2.qasm\n");
        exit(EXIT_FAILURE);
    }

    char fileName1[200];
    char fileName2[200];
    strcpy(fileName1, argv[1]);
    strcpy(fileName2, argv[2]);

    circuit_U = parse_qasm_file(fileName1);
    circuit_V = parse_qasm_file(fileName2);
    optimize_qubit_order(circuit_U, false); // TODO: remove
    optimize_qubit_order(circuit_V, false);
    // TODO: instead: make make qmdd_create_controlled_gate independent of order
    //       of c and t (or make it a wrapper which calls create_multi_cgate)

    // Init
    lace_start(workers, 0);
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, wgt_tab_size, tolerance, wgt_table_type, wgt_norm_strat);

    clock_t cpu_t1 = clock();
    double wall_t1 = wctime();

    circuit_equivalence_check(circuit_U, circuit_V);
        
    clock_t cpu_t2 = clock();
    double wall_t2 = wctime();

    stats.cpu_time  = (double)(cpu_t2 - cpu_t1) / CLOCKS_PER_SEC;
    stats.wall_time = wall_t2 - wall_t1;
    print_stats();

    // Cleanup
    sylvan_quit();
    lace_stop();
    free_quantum_circuit(circuit_U);
    free_quantum_circuit(circuit_V);

    return 0;
}
