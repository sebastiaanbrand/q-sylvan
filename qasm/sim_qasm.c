#include "qsylvan.h"
#include "simple_parser.h"


QMDD apply_gate(QMDD state, quantum_op_t* gate)
{
    // This looks very ugly, but we need some way to match the name of a gate to
    // the GATEID, and since this is C we don't really have easy access to
    // data structures like a dictionary.
    if (strcmp(gate->name, "id") == 0) {
        return state;
    }
    else if (strcmp(gate->name, "x") == 0) {
        return qmdd_gate(state, GATEID_X, gate->target);
    }
    else if (strcmp(gate->name, "y") == 0) {
        return qmdd_gate(state, GATEID_Y, gate->target);
    }
    else if (strcmp(gate->name, "z") == 0) {
        return qmdd_gate(state, GATEID_Z, gate->target);
    }
    else if (strcmp(gate->name, "h") == 0) {
        return qmdd_gate(state, GATEID_H, gate->target);
    }
    else if (strcmp(gate->name, "s") == 0) {
        return qmdd_gate(state, GATEID_S, gate->target);
    }
    else if (strcmp(gate->name, "sdg") == 0) {
        return qmdd_gate(state, GATEID_Sdag, gate->target);
    }
    else if (strcmp(gate->name, "t") == 0) {
        return qmdd_gate(state, GATEID_T, gate->target);
    }
    else if (strcmp(gate->name, "tdg") == 0) {
        return qmdd_gate(state, GATEID_Tdag, gate->target);
    }
    else if (strcmp(gate->name, "sx") == 0) {
        return qmdd_gate(state, GATEID_sqrtX, gate->target);
    }
    else if (strcmp(gate->name, "sxdg") == 0) {
        return qmdd_gate(state, GATEID_sqrtXdag, gate->target);
    }
    else if (strcmp(gate->name, "rx") == 0) {
        return qmdd_gate(state, GATEID_Rx(gate->angle[0]), gate->target);
    }
    else if (strcmp(gate->name, "ry") == 0) {
        return qmdd_gate(state, GATEID_Ry(gate->angle[0]), gate->target);
    }
    else if (strcmp(gate->name, "rz") == 0) {
        return qmdd_gate(state, GATEID_Rz(gate->angle[0]), gate->target);
    }
    else if (strcmp(gate->name, "p") == 0) {
        return qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->target);
    }
    else if (strcmp(gate->name, "u2") == 0) {
        fl_t pi_over_2 = flt_acos(0.0);
        return qmdd_gate(state, GATEID_U(pi_over_2, gate->angle[0], gate->angle[1]), gate->target);
    }
    else if (strcmp(gate->name, "u") == 0) {
        return qmdd_gate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->target);
    }
    else if (strcmp(gate->name, "cx") == 0) {
        return qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "cy") == 0) {
        return qmdd_cgate(state, GATEID_Y, gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "cz") == 0) {
        return qmdd_cgate(state, GATEID_Z, gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "ch") == 0) {
        return qmdd_cgate(state, GATEID_H, gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "csx") == 0) {
        return qmdd_cgate(state, GATEID_sqrtX, gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "crx") == 0) {
        return qmdd_cgate(state, GATEID_Rx(gate->angle[0]), gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "cry") == 0) {
        return qmdd_cgate(state, GATEID_Ry(gate->angle[0]), gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "crz") == 0) {
        return qmdd_cgate(state, GATEID_Rz(gate->angle[0]), gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "cp") == 0) {
        return qmdd_cgate(state, GATEID_Phase(gate->angle[0]), gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "cu") == 0) {
        return qmdd_cgate(state, GATEID_U(gate->angle[0], gate->angle[1], gate->angle[2]), gate->ctrls[0], gate->target);
    }
    else if (strcmp(gate->name, "ccx") == 0) {
        return qmdd_cgate2(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->target);
    }
    else if (strcmp(gate->name, "c3x") == 0) {
        return qmdd_cgate3(state, GATEID_X, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->target);
    }
    else if (strcmp(gate->name, "c3sx") == 0) {
        return qmdd_cgate3(state, GATEID_sqrtX, gate->ctrls[0], gate->ctrls[1], gate->ctrls[2], gate->target);
    }
    else if (strcmp(gate->name, "swap") == 0) {
        // no native SWAP gates in Q-Sylvan
        return qmdd_circuit_swap(state, gate->target, gate->ctrls[0]);
    }
    else if (strcmp(gate->name, "rzz") == 0 ) {
        // no native RZZ gates in Q-Sylvan
        state = qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->target);
        state = qmdd_gate(state, GATEID_Phase(gate->angle[0]), gate->target);
        state = qmdd_cgate(state, GATEID_X, gate->ctrls[0], gate->target);
        return state;
    }
    else {
        fprintf(stderr, "Gate '%s' currently unsupported\n", gate->name);
        return state;
    }
}


QMDD measure(QMDD state, quantum_op_t *measurement)
{
    // TODO
    return state;
}


void simulate_circuit(quantum_circuit_t* circuit)
{
    QMDD state = qmdd_create_all_zero_state(circuit->qreg_size);
    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        if (op->type == op_gate) {
            state = apply_gate(state, op);
        }
        else if (op->type == op_measurement) {
            state = measure(state, op);
        }
        op = op->next;
    }
}


int main(int argc, char *argv[])
{
    quantum_circuit_t* circuit = parse_qasm_file(argv[1]);
    optimize_qubit_order(circuit);
    //print_quantum_circuit(circuit);
    
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_defaults(1LL<<20);

    simulate_circuit(circuit);

    sylvan_quit();
    lace_stop();
    free_quantum_circuit(circuit);

    return 0;
}
