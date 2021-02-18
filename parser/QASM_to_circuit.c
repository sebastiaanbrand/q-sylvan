#include "QASM_to_circuit.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

// TASK_IMPL_3(bool, get_qubits_circuit, char*, token, BDDVAR, n_qubits, BDDVAR*, qubits)
bool get_qubits_circuit(char* token, BDDVAR n_qubits, BDDVAR* qubits)
{
    BDDVAR j = 0;
    strtok(token, "[");
    for(; j < n_qubits; j++) {
        qubits[j] = (BDDVAR)atoi(strtok(NULL, "]"));
        strtok(NULL, "[");
    }
    qubits[j] = (BDDVAR)atoi(strtok(NULL, "]"));
    return true;
}

// TASK_IMPL_2(bool, get_gateid_circuit, char*, gate_str, Gate*, gate_s)
BDDVAR get_gateid_circuit(char* gate_str, Gate* gate_s)
{
    BDDVAR n_qubits = 1;
    while(gate_str[0] == 'c' && strcmp(gate_str, "creg") != 0) {
        gate_str++;
        n_qubits++;
    }
    // Handle gates
    if (strcmp(gate_str, "i") == 0)
        memcpy(gate_s, &gate_I, sizeof(Gate));
    else if (strcmp(gate_str, "h") == 0)
        memcpy(gate_s, &gate_H, sizeof(Gate));
    else if (strcmp(gate_str, "x") == 0)
        memcpy(gate_s, &gate_X, sizeof(Gate));
    else if (strcmp(gate_str, "y") == 0)
        memcpy(gate_s, &gate_Y, sizeof(Gate));
    else if (strcmp(gate_str, "z") == 0)
        memcpy(gate_s, &gate_Z, sizeof(Gate));
    else if (strcmp(gate_str, "s") == 0)
        memcpy(gate_s, &gate_S, sizeof(Gate));
    else if (strcmp(gate_str, "sdg") == 0)
        memcpy(gate_s, &gate_Sd, sizeof(Gate));
    else if (strcmp(gate_str, "t") == 0)
        memcpy(gate_s, &gate_T, sizeof(Gate));
    else if (strcmp(gate_str, "tdg") == 0)
        memcpy(gate_s, &gate_Td, sizeof(Gate));
    else if (strcmp(gate_str, "barrier") == 0)
        memcpy(gate_s, &gate_barrier, sizeof(Gate));
    else if (strcmp(gate_str, "measure") == 0)
        memcpy(gate_s, &gate_measure, sizeof(Gate));
    else {
        gate_str = strtok(gate_str, "(");
        if (strcmp(gate_str, "rx") == 0) {
            gate_str = strtok(NULL, ")");
            memcpy(gate_s, &gate_Rx, sizeof(Gate));
            gate_s->rotation = atof(gate_str);
        } else if (strcmp(gate_str, "ry") == 0) {
            gate_str = strtok(NULL, ")");
            memcpy(gate_s, &gate_Ry, sizeof(Gate));
            gate_s->rotation = atof(gate_str);
        } else if (strcmp(gate_str, "rz") == 0) {
            gate_str = strtok(NULL, ")");
            memcpy(gate_s, &gate_Rz, sizeof(Gate));
            gate_s->rotation = atof(gate_str);
        } else {
            return 0;
        }
    }
    return n_qubits;
}

// TASK_IMPL_5(bool, get_parallel_depth, char*, targets, Gate**, circuit, BDDVAR, n_qubits, BDDVAR*, curr_depth, Gate, gate_s)
bool handle_gate(char* targets, Gate** circuit, BDDVAR n_qubits, BDDVAR* curr_depth, Gate gate_s)
{
    BDDVAR *qubits = malloc((n_qubits+1) * sizeof(BDDVAR));
    if (gate_s.id == gate_barrier.id) {
        *curr_depth += 1;
        for (BDDVAR i = 1; i <= n_qubits; i++)
            circuit[i][*curr_depth] = gate_s;
        BDDVAR *control = malloc((n_qubits) * sizeof(BDDVAR));
        for (BDDVAR i = 0; i < n_qubits; i++) control[i] = i+1;
        gate_s.control = control;
        gate_s.controlSize = n_qubits;
        circuit[0][*curr_depth] = gate_s;
    }
    else {
        get_qubits_circuit(targets, n_qubits, qubits);
        gate_s.controlSize = n_qubits;
        *curr_depth += 1;
        if (n_qubits != 0) {
            BDDVAR *control = malloc((n_qubits) * sizeof(BDDVAR));
            for (BDDVAR i = 0; i < n_qubits; i++) control[i] = qubits[i];
            gate_s.control = control;
            Gate gate_ctrl_s;
            for (BDDVAR i = 0; i < n_qubits; i++) {
                memcpy(&gate_ctrl_s, &gate_ctrl, sizeof(Gate));
                gate_ctrl_s.control = malloc(sizeof(BDDVAR));
                gate_ctrl_s.control[0] = qubits[n_qubits];
                gate_ctrl_s.controlSize = 1;
                circuit[qubits[i]][*curr_depth] = gate_ctrl_s;
            }
        }
        circuit[qubits[n_qubits]][*curr_depth] = gate_s;
    }
    free(qubits);
    return true;
}

// TASK_IMPL_4(bool, handle_line_circuit, char*, line, Gate**, circuit, BDDVAR*, nvars, BDDVAR*, curr_depth)
bool handle_line_circuit(char* line, Gate** circuit, BDDVAR* nvars, BDDVAR* curr_depth)
{
    Gate* gate_s = malloc(sizeof(Gate));
    BDDVAR n_qubits = 0;
    char* gate_str;
    char* targets_str;
    // tokenize string
    gate_str = strtok(line, " ");
    targets_str = strtok(NULL, "");
    // Create all-zero state with "temp" qubits
    if (strstr(gate_str, "qreg") != NULL) {
        get_qubits_circuit(targets_str, 0, nvars);
    }
    else {
        n_qubits = get_gateid_circuit(gate_str, gate_s);
        if (n_qubits != 0) {
            if (gate_s->id == gate_barrier.id) n_qubits = *nvars;
            handle_gate(targets_str, circuit, n_qubits-1, curr_depth, *gate_s);
        }
    }
    free(gate_s);
    return true;
}

Gate** make_circuit(char *filename, bool optimize)
{
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        perror("Error while opening the file.\n");
        return circuit;
    }

    char *line = NULL, *c;
    Gate* realloc_wire;
    circuit = malloc(max_qubits*sizeof(*circuit));
    for (BDDVAR i = 0; i < max_qubits; ++i) {
        circuit[i] = malloc(max_wire * sizeof(Gate));
        for (BDDVAR j = 0; j < max_wire; j++)
            circuit[i][j] = gate_I;
    }
    nvars = malloc(sizeof(BDDVAR));
    curr_depth = malloc(sizeof(BDDVAR));
    *nvars = 0;
    *curr_depth = 0;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, f)) != -1) {
        if (*nvars >= max_qubits)
            perror("Too much qubits, current maximum is 128.");
        while(*curr_depth >= wire_i*max_wire) {
            wire_i++;
            for (BDDVAR i = 0; i < max_qubits; ++i) {
                realloc_wire = realloc(circuit[i], wire_i*max_wire * sizeof *circuit[i]);
                if (realloc_wire == NULL)
                    perror("Memory allocation failed.");
                else
                    circuit[i] = realloc_wire;
                for (BDDVAR j = (wire_i-1)*max_wire; j < wire_i*max_wire; j++)
                    circuit[i][j] = gate_I;
            }
            free(realloc_wire);
        }
        // skip if comment
        if(line[0] == '/' && line[1] == '/')
            continue;
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t') || (*line == '\n'))
            line++;
        // remove empty lines, trailing information after ';'
        c = strchr(line, ';');
        if (c != NULL) {
            line[c - line] = '\0';
            handle_line_circuit(line, circuit, nvars, curr_depth);
        }
    }
    fclose(f);
    if (optimize)
        *curr_depth = optimize_circuit(circuit, *nvars, *curr_depth);
    free(line);
    return circuit;
}

void circuit_exit() {
    for (BDDVAR j = 0; j < *curr_depth; j++) {
        for (BDDVAR i = 0; i < *nvars; i++) {
            if (circuit[i][j].id != gate_I.id) {
                if (circuit[i][j].control != NULL)
                    free(circuit[i][j].control);
            }
        }
    }
    for (BDDVAR i = 0; i < max_qubits; i++)
        free(circuit[i]);
    free(circuit);
    free(nvars);
    free(curr_depth);
}

BDDVAR optimize_circuit(Gate** circuit, BDDVAR nvars, BDDVAR depth)
{
    BDDVAR depth2;
    for (BDDVAR q = 0; q < nvars; q++) {
        for (BDDVAR depth1 = 0; depth1 < depth; depth1++) {
            if (circuit[q][depth1].id != gate_I.id && circuit[q][depth1].id != gate_ctrl.id) {
                depth2 = depth1+1;
                while (circuit[q][depth2].id == gate_I.id || circuit[q][depth2].id == gate_barrier.id) {
                    depth2++;
                }
                if (find_palindromes(circuit, q, depth1, depth2))
                    optimize_circuit_p(circuit, q, depth1, depth2, nvars, depth);
            }
        }
    }
    return reduce_circuit(circuit, nvars, depth);
}

void optimize_circuit_p(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2, BDDVAR nvars, BDDVAR depth)
{
    bool found = true;
    while (found) {
        found = false;
        remove_gates(circuit, q, depth1, depth2);
        do {
            if (depth1 == 0) return;
            depth1--;
        }
        while (circuit[q][depth1].id == gate_barrier.id);
        do {
            if (depth2 == depth) return;
            depth2++;
        }
        while (circuit[q][depth2].id == gate_barrier.id);
        
        for (BDDVAR q = 0; q < nvars; q++) {
            if (circuit[q][depth1].id != gate_I.id && circuit[q][depth1].id != gate_ctrl.id) {
                if(find_palindromes(circuit, q, depth1, depth2))
                    optimize_circuit_p(circuit, q, depth1, depth2, nvars, depth);
            }
        }
    }
}

bool find_palindromes(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    if (circuit[q][depth1].id != circuit[q][depth2].id || fmod(circuit[q][depth1].rotation + circuit[q][depth2].rotation, 1.f) != 0.f)
        return false;
    if (circuit[q][depth1].controlSize == 0)
        return true;
    else if(circuit[q][depth1].controlSize != circuit[q][depth2].controlSize)
        return false;

    for (BDDVAR i = 0; i < circuit[q][depth1].controlSize; i++) {
        if (circuit[q][depth1].control[i] != circuit[q][depth2].control[i])
            return false;
    }
    for (BDDVAR i = 0; i < circuit[q][depth1].controlSize; i++) {
        for(BDDVAR j = depth1+1; j < depth2; j++) {
            if(circuit[circuit[q][depth1].control[i]][j].id != gate_I.id)
                return false;
        }
    }
    return true;
}

void remove_gates(Gate** circuit, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    for (BDDVAR i = 0; i < circuit[q][depth1].controlSize; i++) {
        circuit[circuit[q][depth1].control[i]][depth1] = gate_I;
        circuit[circuit[q][depth1].control[i]][depth2] = gate_I;
    }
    circuit[q][depth1] = gate_I;
    circuit[q][depth2] = gate_I;
}

BDDVAR reduce_circuit(Gate** circuit, BDDVAR nvars, BDDVAR depth)
{
    BDDVAR k;
    // Reduce depth by moving gates to the left if possible
    for (BDDVAR j = 1; j <= depth; j++) {
        for (BDDVAR i = 0; i < nvars; i++) {
            if (circuit[i][j].id != gate_I.id && circuit[i][j].id != gate_ctrl.id)
                reduce_gate(circuit, i, j);
        }
    }
    k = 0;
    // Get new depth
    for (BDDVAR i = 0; i < nvars; i++) {
        for (BDDVAR j = depth; j-- > 0;) {
            if (circuit[i][j].id != gate_I.id) {
                k = j < k ? k : j;
                break;
            }
        }
    }
    return k;
}

void reduce_gate(Gate** circuit, BDDVAR target, BDDVAR depth)
{
    BDDVAR curr;
    BDDVAR reduce = get_reduce_depth(circuit, target, depth);
    for (BDDVAR i = 0; i < circuit[target][depth].controlSize; i++) {
        curr = get_reduce_depth(circuit, circuit[target][depth].control[i], depth);
        if (curr > reduce) reduce = curr;
    }
    if (reduce - depth <= 0) return;
    for (BDDVAR i = 0; i < circuit[target][depth].controlSize; i++) {
        circuit[circuit[target][depth].control[i]][reduce] = circuit[circuit[target][depth].control[i]][depth];
        circuit[circuit[target][depth].control[i]][depth] = gate_I;
    }
    circuit[target][reduce] = circuit[target][depth];
    circuit[target][depth] = gate_I;
}

BDDVAR get_reduce_depth(Gate** circuit, BDDVAR target, BDDVAR depth)
{
    while (depth > 0 && (circuit[target][depth-1].id == gate_I.id))
        depth--;
    return depth;
}

void print_circuit(Gate** circuit, BDDVAR* nvars, BDDVAR* curr_depth)
{
    for (BDDVAR i = 0; i < *nvars; i++) {
        for (BDDVAR j = 0; j < *curr_depth+1; j++) {
            printf("-%s",circuit[i][j].gateSymbol);
        }
        printf("-\n");
    }
}

int main(int argc, char *argv[])
{
    // check for file
    char *filename = "";
    bool optimize = false;
    int opt;
    while((opt = getopt(argc, argv, "f:o")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 'o':
                optimize = true;
                break;
        }
    }
    if(strcmp(filename, "") == 0)
    {
        printf("Give filename of qasm file.\n");
        return 1;
    }

    make_circuit(filename, optimize);
    print_circuit(circuit, nvars, curr_depth);
    circuit_exit();
    return 0;
}