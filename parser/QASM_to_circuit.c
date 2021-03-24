#include "QASM_to_circuit.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

C_struct make_c_struct(char *filename, bool optimize)
{
    // Open the QASM file (or give error)
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        perror("File error");
        exit(0);
    }

    C_struct c_s = c_struct_default;
    // Allocate space for the circuit
    c_s.circuit = malloc(c_s.max_qubits*sizeof(c_s.circuit));
    for (BDDVAR i = 0; i < c_s.max_qubits; ++i) {
        c_s.circuit[i] = malloc(c_s.max_wire * sizeof(Gate));
        for (BDDVAR j = 0; j < c_s.max_wire; j++)
            c_s.circuit[i][j] = gate_I;
    }

    char* line = NULL;
    char* c;
    size_t len = 0;
    ssize_t read;
    // reallocate_wire(&c_s);

    while ((read = getline(&line, &len, f)) != -1) {
        if (c_s.nvars >= c_s.max_qubits) {
            perror("Too much qubits, current maximum is 128.");
            exit(0);
        }
        if(c_s.depth == c_s.max_wire-1)
            reallocate_wire(&c_s);
        // skip if comment
        if(line[0] == '/' && line[1] == '/')
            continue;
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t'))
            line++;
        // remove empty lines, trailing information after ';'
        c = strchr(line, ';');
        if (c != NULL) {
            line[c - line] = '\0';
            handle_line_c_struct(line, &c_s);
        }
    }
    free(line);
    fclose(f);

    if (optimize)
        optimize_c_struct(&c_s);
    reduce_c_struct(&c_s);
    return c_s;
}

void reallocate_wire(C_struct* c_s)
{
    BDDVAR incr = 1024;
    c_s->max_wire += incr;
    for (BDDVAR i = 0; i < c_s->max_qubits; i++) {
        c_s->circuit[i] = realloc(c_s->circuit[i], c_s->max_wire * sizeof(Gate));
        if (c_s->circuit[i] == NULL) {
            perror("Memory allocation failed.");
            delete_c_struct(c_s);
            exit(0);
        }
        for (BDDVAR j = c_s->max_wire-incr; j < c_s->max_wire; j++)
            c_s->circuit[i][j] = gate_I;
    }
}

void delete_c_struct(C_struct* c_s)
{
    for (BDDVAR j = 0; j < c_s->max_wire; j++) {
        for (BDDVAR i = 0; i < c_s->max_qubits; i++) {
            if (c_s->circuit[i][j].id != gate_I.id) {
                if (c_s->circuit[i][j].control != NULL || c_s->circuit[i][j].controlSize != 0)
                    free(c_s->circuit[i][j].control);
            }
        }
    }
    for (BDDVAR i = 0; i < c_s->max_qubits; i++)
        free(c_s->circuit[i]);
    free(c_s->circuit);
}

bool handle_line_c_struct(char* line, C_struct* c_s)
{
    Gate gate_s;
    BDDVAR n_qubits = 0;
    // tokenize string
    char* gate_str = strtok(line, " ");
    char* targets_str = strtok(NULL, "");

    // Create all-zero state with "temp" qubits
    if (strstr(gate_str, "qreg") != NULL)
        get_qubits_c_struct(targets_str, 0, &c_s->nvars);
    else {
        n_qubits = get_gateid_c_struct(gate_str, &gate_s);
        if (n_qubits != 0) {
            if (gate_s.id == gate_barrier.id)
               handle_barrier(c_s, gate_s);
            else {
                BDDVAR *qubits = malloc((n_qubits) * sizeof(BDDVAR));
                get_qubits_c_struct(targets_str, n_qubits-1, qubits);
                handle_gate(c_s, n_qubits-1, qubits, gate_s);
                free(qubits);
            }
        }
    }
    return true;
}

bool get_qubits_c_struct(char* token, BDDVAR n_qubits, BDDVAR* qubits)
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

BDDVAR get_gateid_c_struct(char* gate_str, Gate* gate_s)
{
    BDDVAR n_qubits = 1;
    while(gate_str[0] == 'c' && strcmp(gate_str, "creg") != 0) {
        gate_str++;
        n_qubits++;
    }
    // Handle gates
    if (strcmp(gate_str, "i") == 0)
        *gate_s = gate_I;
    else if (strcmp(gate_str, "x") == 0)
        *gate_s = gate_X;
    else if (strcmp(gate_str, "y") == 0)
        *gate_s = gate_Y;
    else if (strcmp(gate_str, "z") == 0)
        *gate_s = gate_Z;
    else if (strcmp(gate_str, "h") == 0)
        *gate_s = gate_H;
    else if (strcmp(gate_str, "sx") == 0)
        *gate_s = gate_sX;
    else if (strcmp(gate_str, "sy") == 0)
        *gate_s = gate_sY;
    else if (strcmp(gate_str, "s") == 0)
        *gate_s = gate_S;
    else if (strcmp(gate_str, "sdg") == 0)
        *gate_s = gate_Sd;
    else if (strcmp(gate_str, "t") == 0)
        *gate_s = gate_T;
    else if (strcmp(gate_str, "tdg") == 0)
        *gate_s = gate_Td;
    else if (strcmp(gate_str, "barrier") == 0)
        *gate_s = gate_barrier;
    else if (strcmp(gate_str, "measure") == 0)
        *gate_s = gate_measure;
    else {
        gate_str = strtok(gate_str, "(");
        if (strcmp(gate_str, "rx") == 0) {
            gate_str = strtok(NULL, ")");
            *gate_s = gate_Rx;
            gate_s->rotation = atof(gate_str);
        } else if (strcmp(gate_str, "ry") == 0) {
            gate_str = strtok(NULL, ")");
            *gate_s = gate_Ry;
            gate_s->rotation = atof(gate_str);
        } else if (strcmp(gate_str, "rz") == 0) {
            gate_str = strtok(NULL, ")");
            *gate_s = gate_Rz;
            gate_s->rotation = atof(gate_str);
        } else {
            return 0;
        }
    }
    return n_qubits;
}

void copy_Gate(Gate src, Gate* dst)
{
    memcpy(dst, &src, sizeof(Gate));
    *dst->gateSymbol = *src.gateSymbol;
}

void handle_barrier(C_struct* c_s, Gate gate_s)
{
    c_s->depth++;
    for (BDDVAR i = 0; i < c_s->nvars-1; i++)
        c_s->circuit[i][c_s->depth] = gate_s;
    BDDVAR *control = malloc((c_s->nvars-1) * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s->nvars-1; i++) control[i] = i;
    gate_s.control = control;
    gate_s.controlSize = c_s->nvars-1;
    c_s->circuit[c_s->nvars-1][c_s->depth] = gate_s;
}

void handle_gate(C_struct* c_s, BDDVAR n_qubits, BDDVAR* qubits, Gate gate_s)
{
    gate_s.controlSize = n_qubits;
    c_s->depth += 1;
    if (n_qubits != 0) {
        BDDVAR *control = malloc((n_qubits) * sizeof(BDDVAR));
        for (BDDVAR i = 0; i < n_qubits; i++) control[i] = qubits[i];
        gate_s.control = control;
        Gate gate_ctrl_s;
        bool controlled = false;
        BDDVAR j = 0;
        for (BDDVAR i = 0; i < c_s->nvars; i++) {
            if (i == qubits[n_qubits]) break;
            if (qubits[j] == i) {
                controlled = true;
                j++;
                gate_ctrl_s = gate_ctrl;
                gate_ctrl_s.control = malloc(sizeof(BDDVAR));
                gate_ctrl_s.control[0] = qubits[n_qubits];
                gate_ctrl_s.controlSize = 1;
                c_s->circuit[i][c_s->depth] = gate_ctrl_s;
            }
            else if (controlled)
                c_s->circuit[i][c_s->depth] = gate_ctrl_c;
        }
    }
    c_s->circuit[qubits[n_qubits]][c_s->depth] = gate_s;
}

void optimize_c_struct(C_struct* c_s)
{
    BDDVAR depth2;
    for (BDDVAR q = 0; q < c_s->nvars; q++) {
        for (BDDVAR depth1 = 0; depth1 < c_s->depth; depth1++) {
            if (c_s->circuit[q][depth1].id != gate_I.id && c_s->circuit[q][depth1].id != gate_ctrl.id) {
                depth2 = depth1+1;
                while (c_s->circuit[q][depth2].id == gate_I.id || c_s->circuit[q][depth2].id == gate_barrier.id) {
                    if (depth2 >= c_s->depth) break;
                    depth2++;
                }
                if (find_palindromes(c_s, q, depth1, depth2))
                    optimize_c_struct_p(c_s, q, depth1, depth2);
            }
        }
    }
}

void optimize_c_struct_p(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    bool found = true;
    while (found) {
        found = false;
        remove_gates(c_s, q, depth1, depth2);
        do {
            if (depth1 == 0) return;
            depth1--;
        }
        while (c_s->circuit[q][depth1].id == gate_barrier.id || c_s->circuit[q][depth1].id == gate_I.id);

        do {
            if (depth2 >= c_s->depth) return;
            depth2++;
        }
        while (c_s->circuit[q][depth2].id == gate_barrier.id || c_s->circuit[q][depth2].id == gate_I.id);

        for (BDDVAR q = 0; q < c_s->nvars; q++) {
            if (c_s->circuit[q][depth1].id != gate_I.id && c_s->circuit[q][depth1].id != gate_ctrl.id) {
                if(find_palindromes(c_s, q, depth1, depth2))
                    optimize_c_struct_p(c_s, q, depth1, depth2);
            }
        }
    }
}

bool find_palindromes(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    if (c_s->circuit[q][depth1].id != c_s->circuit[q][depth2].id || 
            c_s->circuit[q][depth1].rotation != c_s->circuit[q][depth2].rotation || 
            c_s->circuit[q][depth1].rotation != -c_s->circuit[q][depth2].rotation)
        return false;
    if (c_s->circuit[q][depth1].controlSize == 0)
        return true;
    else if(c_s->circuit[q][depth1].controlSize != c_s->circuit[q][depth2].controlSize)
        return false;

    for (BDDVAR i = 0; i < c_s->circuit[q][depth1].controlSize; i++) {
        if (c_s->circuit[q][depth1].control[i] != c_s->circuit[q][depth2].control[i])
            return false;
    }
    for (BDDVAR i = 0; i < c_s->circuit[q][depth1].controlSize; i++) {
        for(BDDVAR j = depth1+1; j < depth2; j++) {
            if(c_s->circuit[c_s->circuit[q][depth1].control[i]][j].id != gate_I.id)
                return false;
        }
    }
    return true;
}

void remove_gates(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    for (BDDVAR i = 0; i < c_s->circuit[q][depth1].controlSize; i++) {
        c_s->circuit[c_s->circuit[q][depth1].control[i]][depth1] = gate_I;
        c_s->circuit[c_s->circuit[q][depth1].control[i]][depth2] = gate_I;
    }
    c_s->circuit[q][depth1] = gate_I;
    c_s->circuit[q][depth2] = gate_I;
}

void reduce_c_struct(C_struct* c_s)
{
    BDDVAR k;
    // Reduce depth by moving gates to the left if possible
    for (BDDVAR j = 1; j <= c_s->depth; j++) {
        for (BDDVAR i = 0; i < c_s->nvars; i++) {
            if ((c_s->circuit[i][j].id != gate_I.id && 
                c_s->circuit[i][j].id != gate_ctrl.id && 
                c_s->circuit[i][j].id != gate_ctrl_c.id &&
                c_s->circuit[i][j].id != gate_barrier.id) || 
                (c_s->circuit[i][j].id == gate_barrier.id && i == c_s->nvars-1))
                reduce_gate(c_s, i, j);
        }
    }
    k = 0;
    // Get new depth
    for (BDDVAR i = 0; i < c_s->nvars; i++) {
        for (BDDVAR j = c_s->depth; j-- > 0;) {
            if (c_s->circuit[i][j].id != gate_I.id) {
                k = j < k ? k : j;
                break;
            }
        }
    }
    c_s->depth = k+1;
}

void reduce_gate(C_struct* c_s, BDDVAR target, BDDVAR depth)
{
    BDDVAR curr;
    BDDVAR reduce = get_reduce_depth(c_s, target, depth);
    BDDVAR *control = c_s->circuit[target][depth].control;
    if (control != NULL) {
        for (BDDVAR i = control[0]; i < target; i++) {
            curr = get_reduce_depth(c_s, i, depth);
            if (curr > reduce) reduce = curr;
        }
        if (reduce - depth <= 0) return;
        for (BDDVAR i = control[0]; i < target; i++) {
            c_s->circuit[i][reduce] = c_s->circuit[i][depth];
            c_s->circuit[i][depth] = gate_I;
        }
    }
    if (reduce - depth <= 0) return;
    c_s->circuit[target][reduce] = c_s->circuit[target][depth];
    c_s->circuit[target][depth] = gate_I;
}

BDDVAR get_reduce_depth(C_struct* c_s, BDDVAR target, BDDVAR depth)
{
    while (depth > 0 && (c_s->circuit[target][depth-1].id == gate_I.id))
        depth--;
    return depth;
}

void print_c_struct(C_struct c_s, bool vertical, bool show_rotation)
{
    bool has_rotation = false;
    bool negative_rotation = false;
    if (vertical) {
        for (BDDVAR j = 0; j < c_s.depth; j++) {
            for (BDDVAR i = 0; i < c_s.nvars; i++)
                printf(" |  ");
            printf("\n");
            for (BDDVAR i = 0; i < c_s.nvars; i++) {
                if (c_s.circuit[i][j].id == gate_barrier.id)
                    printf("----");
                else {
                    if (c_s.circuit[i][j].gateSymbol[0] == '-')
                        printf(" |");
                    else
                        printf(" %c",c_s.circuit[i][j].gateSymbol[0]);
                    if (c_s.circuit[i][j].gateSymbol[1] != '-')
                        printf("%c ",c_s.circuit[i][j].gateSymbol[1]);
                    else
                        printf("  ");
                }
            }
            printf("\n");
        }
    }
    else {
        for (BDDVAR i = 0; i < c_s.nvars; i++) {
            for (BDDVAR j = 0; j < c_s.depth; j++) {
                has_rotation = false;
                negative_rotation = false;
                for (BDDVAR k = 0; k < c_s.nvars; k++) {
                    if(c_s.circuit[k][j].id == 11 || c_s.circuit[k][j].id == 12 || c_s.circuit[k][j].id == 13) has_rotation = true;
                    if(c_s.circuit[k][j].rotation < 0) negative_rotation = true;
                }
                if (has_rotation && show_rotation) {
                    if(c_s.circuit[i][j].id == 11 || c_s.circuit[i][j].id == 12 || c_s.circuit[i][j].id == 13)
                        printf("-%s(%.4lf)",c_s.circuit[i][j].gateSymbol,roundf(c_s.circuit[i][j].rotation*10000)/10000);
                    else {
                        if (negative_rotation)
                            printf("-%s---------",c_s.circuit[i][j].gateSymbol);
                        else
                            printf("-%s--------",c_s.circuit[i][j].gateSymbol);
                    }
                }
                else
                    printf("-%s",c_s.circuit[i][j].gateSymbol);
            }
            printf("-\n");
        }
    }
}