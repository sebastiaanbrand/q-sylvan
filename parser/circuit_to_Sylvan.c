#include "circuit_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"
#include <time.h> 

#include <sys/time.h>
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

// TASK_IMPL_4(bool, measure, QDD, qdd, bool*, measurements, BDDVAR, nvars, BDDVAR, shots)
bool measure(QDD qdd, bool* measurements, BDDVAR nvars, BDDVAR shots, bool show)
{
    // Run the circuit x times, where x == shots
    bool *ms = malloc(nvars * sizeof(bool));
    double *p = malloc(nvars * sizeof(double));
    int *hits = malloc(pow(2,nvars) * sizeof(int));
    for(int i = 0; i < pow(2,nvars); i++) { hits[i] = 0; }
    for(BDDVAR i = 0; i < shots; i++) { 
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

    // Print reformatted circuit results
    if (show) {
        for(int k = 0; k < pow(2, sum); k++) {
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
    free(ms);
    free(p);
    free(hits);
    free(probs);
    return true;
}

TASK_IMPL_1(BDDVAR, get_gate_id, Gate, gate)
// BDDVAR get_gate_id(Gate gate)
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

TASK_IMPL_3(QDD, apply_gate, QDD, qdd, Gate, gate, BDDVAR, i)
{
    BDDVAR gate_id = get_gate_id(gate);
    if (gate.controlSize == 0)
        qdd = qdd_gate(qdd, gate_id, i);
    else if (gate.controlSize == 1)
        qdd = qdd_cgate(qdd, gate_id, gate.control[0], i);
    else if (gate.controlSize == 2)
        qdd = qdd_cgate2(qdd, gate_id, gate.control[0], gate.control[1], i);
    else if (gate.controlSize == 3)
        qdd = qdd_cgate3(qdd, gate_id, gate.control[0], gate.control[1], gate.control[2], i);
    return qdd;
}

bool run_circuit_matrix(C_struct c_s, BDDVAR shots, bool show)
{
    LACE_ME;
    QDD vec = qdd_create_all_zero_state(c_s.nvars);
    QDD qdd, ctrl_qdd, col_qdd, new_qdd;
    Gate gate = gate_barrier;
    uint32_t gate_id;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    int* options = malloc(c_s.nvars * sizeof(int));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    uint32_t* column_of_gates = malloc(c_s.nvars*sizeof(uint32_t));
    qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
    for (BDDVAR j = 0; j < c_s.depth; j++) {
        new_qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            column_of_gates[i] = GATEID_I;
            gate = c_s.circuit[i][j];
            if (gate.id == gate_barrier.id) {
                vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
                qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
                break;
            }
            else if (gate.id == gate_measure.id)
                measurements[i] = true;
            else if (gate.id != gate_ctrl.id && gate.id != gate_I.id) {
                gate_id = get_gate_id(c_s.circuit[i][j]);
                if (gate.controlSize == 0)
                    column_of_gates[i] = gate_id;
                else {
                    column_of_gates[i] = GATEID_I;
                    for (BDDVAR k = 0; k < c_s.nvars; k++) options[k] = -1;
                    for (BDDVAR k = 0; k < gate.controlSize; k++) options[gate.control[k]] = 1;
                    options[i] = 2;
                    ctrl_qdd = qdd_create_multi_cgate_rec(c_s.nvars, options, gate_id, 0);
                    new_qdd = qdd_matmat_mult(ctrl_qdd, new_qdd, c_s.nvars);
                }
            }
        }
        if (gate.id != gate_barrier.id) {
            col_qdd = qdd_create_single_qubit_gates(c_s.nvars, column_of_gates);
            new_qdd = qdd_matmat_mult(col_qdd, new_qdd, c_s.nvars);
            qdd = qdd_matmat_mult(new_qdd, qdd, c_s.nvars);
        }
        // printf("%d/%d: %ld\n", (int)floor(4*pow(c_s.nvars,2)/5),(int)pow(c_s.nvars,2),qdd_countnodes(qdd));
        if (qdd_countnodes(qdd) > floor(4*pow(c_s.nvars,2)/5)) {
            // printf("%d-",j);
            vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
            qdd = qdd_create_single_qubit_gates_same(c_s.nvars, GATEID_I);
        }
    }
    printf("%ld %ld\n", qdd_countnodes(qdd),qdd_countnodes(vec));
    vec = qdd_matvec_mult(qdd, vec, c_s.nvars);
    measure(vec, measurements, c_s.nvars, shots, show);
    return true;
}

bool run_c_struct(C_struct c_s, BDDVAR shots, bool show)
{
    LACE_ME;

    Gate gate;
    bool* measurements = malloc(c_s.nvars * sizeof(bool));
    for (BDDVAR i = 0; i < c_s.nvars; i++) measurements[i] = false;
    QDD qdd = qdd_create_all_zero_state(c_s.nvars);

    for (BDDVAR j = 0; j < c_s.depth; j++) {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            gate = c_s.circuit[i][j];
            if (gate.id == gate_barrier.id || gate.id == gate_ctrl.id || gate.id == gate_I.id)
                continue;
            else if (gate.id == gate_measure.id) {
                measurements[i] = true;
                continue;
            }
            else
                qdd = apply_gate(qdd, gate, i);
        }
    }
    measure(qdd, measurements, c_s.nvars, shots, show);
    free(measurements);
    return true;
}

int main(int argc, char *argv[])
{
    // LACE_ME;
    // check for file
    char *filename = "../plotting/test_gates.txt";
    BDDVAR shots = 100;
    BDDVAR seed = 100;
    bool matrix = false;
    int opt;
    while((opt = getopt(argc, argv, "f:s:r:m")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 's':
                shots = atoi(optarg);
                break;
            case 'r':
                seed = atoi(optarg);
                break;
            case 'm':
                matrix = true;
                break;
        }
    }
    srand(seed);
    if(strcmp(filename, "") == 0)
    {
        printf("Give filename of qasm file.\n");
        return 1;
    }
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<16, -1, true);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    C_struct c_s = make_c_struct(filename, true);
    double begin = wctime();
    if (matrix)
        run_circuit_matrix(c_s, shots, false);
    else
        run_c_struct(c_s, shots, false);
    double end = wctime();
    double time_spent = end - begin;
    printf("time: %f\n", (time_spent));
    delete_c_struct(&c_s);
    sylvan_quit();
    lace_exit();

    return 0;
}