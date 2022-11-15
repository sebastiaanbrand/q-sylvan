#include "QASM_to_circuit.h"

C_struct make_c_struct(char *filename, bool optimize)
{
    // Open the QASM file (or give error)
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Unable to open QASM file.\n");
        exit(EXIT_FAILURE);
    }
    // Create a default circuit struct
    C_struct c_s = c_struct_default;
    // Allocate space for the circuit
    c_s.circuit = malloc(c_s.max_qubits*sizeof(c_s.circuit));
    for (BDDVAR i = 0; i < c_s.max_qubits; ++i) {
        c_s.circuit[i] = malloc(c_s.max_wire * sizeof(Gate));
        for (BDDVAR j = 0; j < c_s.max_wire; j++)
            c_s.circuit[i][j] = gate_I;
    }
    // Initialise variables
    char* line = NULL;
    char* c;
    size_t len = 0;
    ssize_t read;
    int index = 0;
    // Read all lines from the file
    while ((read = getline(&line, &len, f)) != -1) {
        index++;
        // If the circuit takes more than 128 qubits, give error
        if (c_s.qubits >= c_s.max_qubits) {
            fprintf(stderr, "Too many qubits, current maximum is 128.\n");
            exit(EXIT_FAILURE);
        }
        // If the current maximum depth is reached, increase the depth and reallocate
        if(c_s.depth == c_s.max_wire-1)
            reallocate_wire(&c_s);
        // skip if comment
        if(line[0] == '/' && line[1] == '/')
            continue;
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t'))
            line++;
        // remove empty lines and trailing information after ';'
        c = strchr(line, ';');
        if (c != NULL) {
            line[c - line] = '\0';
            // Handle line
            if(!handle_line_c_struct(line, &c_s)) {
                fprintf(stderr, "Line %d: invalid QASM code.\n", index);
                exit(EXIT_FAILURE);
            }
        }
    }
    // Free variables
    free(line);
    fclose(f);

    // Optimize if needed
    if (optimize)
        optimize_c_struct(&c_s);
    // Reduce circuit depth by parallelising
    reduce_c_struct(&c_s);
    // Reduce circuit memory by copying the circuit struct (until depth)
    C_struct final_c_s = copy_c_struct(&c_s);
    return final_c_s;
}

void reallocate_wire(C_struct* c_s)
{
    // Increase maximum depth
    BDDVAR incr = 1024;
    c_s->max_wire += incr;
    // Reallocate each wire
    for (BDDVAR i = 0; i < c_s->max_qubits; i++) {
        c_s->circuit[i] = realloc(c_s->circuit[i], c_s->max_wire * sizeof(Gate));
        // If reallocation failed, throw error
        if (c_s->circuit[i] == NULL) {
            fprintf(stderr, "Memory allocation failed.");
            delete_c_struct(c_s);
            exit(EXIT_FAILURE);
        }
        // Fill new space with identity gates
        for (BDDVAR j = c_s->max_wire-incr; j < c_s->max_wire; j++)
            c_s->circuit[i][j] = gate_I;
    }
}

C_struct copy_c_struct(C_struct* c_s)
{
    // Create a default circuit struct
    C_struct new_c_s = c_struct_default;
    // Allocate space for the circuit
    new_c_s.circuit = malloc(c_s->qubits*sizeof(c_s->circuit));
    for (BDDVAR i = 0; i < c_s->qubits; ++i) {
        new_c_s.circuit[i] = malloc(c_s->depth * sizeof(Gate));
        for (BDDVAR j = 0; j < c_s->depth; j++)
            new_c_s.circuit[i][j] = c_s->circuit[i][j];
    }
    new_c_s.qubits = c_s->qubits;
    new_c_s.bits = c_s->bits;
    new_c_s.depth = c_s->depth;
    new_c_s.max_qubits = c_s->qubits;
    new_c_s.max_wire = c_s->depth;
    return new_c_s;
}

void delete_c_struct(C_struct* c_s)
{
    // Loop over all positions in the circuit
    for (BDDVAR j = 0; j < c_s->max_wire; j++) {
        for (BDDVAR i = 0; i < c_s->max_qubits; i++) {
            // If the gate has control qubits, free the list of indices
            if (c_s->circuit[i][j].id != gate_I.id) {
                if (c_s->circuit[i][j].control != NULL || c_s->circuit[i][j].controlSize != 0)
                    free(c_s->circuit[i][j].control);
            }
        }
    }
    // Free each wire
    for (BDDVAR i = 0; i < c_s->max_qubits; i++)
        free(c_s->circuit[i]);
    // Free the remaining data from the circuit struct
    free(c_s->circuit);
}

bool handle_line_c_struct(char* line, C_struct* c_s)
{
    // Initialise variables
    Gate gate_s;
    int n_qubits = 0;
    // tokenize string
    char* gate_str = strtok(line, " ");
    char* targets_str = strtok(NULL, "");
    // Create all-zero state with given size
    if (strstr(gate_str, "qreg") != NULL) {
        if(!get_qubits_c_struct(targets_str, 0, &c_s->qubits))
            return false;
    }
    else if (strstr(gate_str, "creg") != NULL) {
        if(!get_qubits_c_struct(targets_str, 0, &c_s->bits))
            return false;
    }
    else {
        // Get the gate struct corresponding to the gate command and the qubit count
        n_qubits = get_gateid_c_struct(gate_str, &gate_s);
        // If -1 is returned, an if statement is found
        if (n_qubits == -1) {
            char* reg = gate_str;
            gate_str = strtok(targets_str, " ");
            targets_str = strtok(NULL, "");
            // Get the gate struct corresponding to the gate command and the qubit count
            n_qubits = get_gateid_c_struct(gate_str, &gate_s);
            reg = strtok(reg, "=");
            BDDVAR expected = 1;
            gate_str = strtok(NULL, "=");
            if (gate_str != NULL) {
                gate_str = strtok(gate_str, ")");
                expected = (BDDVAR)atoi(gate_str);
            }
            reg = strtok(reg, "[");
            reg = strtok(NULL, "]");
            // If we can find a bracket, the if statement is on a single bit
            if (reg != NULL)
                gate_s.classical_control = (int)atoi(reg);
            gate_s.classical_expect = expected;
        }
        // If n_qubits is 0 there is an unknown command in <line> and if it is -1 the if statement failed
        if (n_qubits > 0) {
            // Handle barrier
            if (gate_s.id == gate_barrier.id)
               handle_barrier(c_s, gate_s);
            else {
                // If it is a measure gate, at a control to store the classical target
                if (gate_s.id == gate_measure.id)
                    n_qubits++;
                // Get the target qubit(s) and handle the gate
                BDDVAR *qubits = malloc((n_qubits) * sizeof(BDDVAR));
                if(!get_qubits_c_struct(targets_str, n_qubits-1, qubits))
                    return false;
                // Check if all qubits are within range of the circuit
                for (int i = 0; i < n_qubits; i++) {
                    if(qubits[i] >= c_s->qubits)
                        return false;
                }
                // Handle measure separately
                if (gate_s.id == gate_measure.id)
                    handle_measure(c_s, gate_s, qubits);
                // Store the gate and free variables
                else
                    handle_gate(c_s, n_qubits-1, qubits, gate_s);
                free(qubits);
            }
        }
    }
    return true;
}

bool get_qubits_c_struct(char* token, BDDVAR n_qubits, BDDVAR* qubits)
{
    // Initialise variables
    BDDVAR j = 0;
    char* index;
    // Find first opening bracket
    strtok(token, "[");
    for(; j < n_qubits; j++) {
        index = strtok(NULL, "]");
        // If NULL, token does not contain enough indices
        if (index == NULL)
            return false;
        // Store index between brackets
        qubits[j] = (BDDVAR)atoi(index);
        // Get next opening bracket
        strtok(NULL, "[");
    }
    index = strtok(NULL, "]");
    // If NULL, token does not contain enough indices
    if (index == NULL)
        return false;
    // Store last index between brackets
    qubits[j] = (BDDVAR)atoi(index);
    return true;
}

int get_gateid_c_struct(char* gate_str, Gate* gate_s)
{
    // Count the controls, starting by 1 (there is always a target)
    // Skip creg, not needed (TODO: maybe implement later?)
    int n_qubits = 1;
    while(gate_str[0] == 'c' && strcmp(gate_str, "creg") != 0) {
        gate_str++;
        n_qubits++;
    }
    // Handle gates
    if (gate_str[0] == 'i' && gate_str[1] == 'f')
        return -1;
    else if (strcmp(gate_str, "i") == 0)
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
        // Handle special rotation gates
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
            // Unknown gate
            if (strcmp(gate_str, "OPENQASM") == 0) return 0;
            if (strcmp(gate_str, "include") == 0) return 0;
            fprintf(stderr, "Error: gate %s currently not supported\n", gate_str);
            exit(1);
        }
    }
    return n_qubits;
}

void handle_measure(C_struct* c_s, Gate gate_s, BDDVAR* qubits)
{
    // Increment the depth since were adding a column
    c_s->depth++;
    // Set target register as control
    gate_s.controlSize = 1;
    gate_s.control = malloc(sizeof(BDDVAR));
    gate_s.control[0] = qubits[1];
    // Set measure on correct qubit
    c_s->circuit[qubits[0]][c_s->depth] = gate_s;
}

void handle_barrier(C_struct* c_s, Gate gate_s)
{
    // Increment the depth since were adding a column
    c_s->depth++;
    // Place barriers on all qubits except the last
    for (BDDVAR i = 0; i < c_s->qubits-1; i++)
        c_s->circuit[i][c_s->depth] = gate_s;
    // Create a list of control indices and store in barrier gate
    // These controls are needed for removing consecutive barriers
    BDDVAR *control = malloc((c_s->qubits-1) * sizeof(BDDVAR));
    for (BDDVAR i = 0; i < c_s->qubits-1; i++) control[i] = i;
    gate_s.control = control;
    gate_s.controlSize = c_s->qubits-1;
    // Set controlled barrier gate as last gate
    c_s->circuit[c_s->qubits-1][c_s->depth] = gate_s;
}

void swap(BDDVAR* a, BDDVAR* b)
{
    BDDVAR t = *a;
    *a = *b;
    *b = t;
}

BDDVAR partition(BDDVAR* qubits, int low, int high)
{
    BDDVAR pivot = qubits[high];
    int i = (low - 1);
    for (int j = low; j < high; j++) {
        if (qubits[j] < pivot) {
            i++;
            swap(&qubits[i], &qubits[j]);
        }
    }
    i++;
    swap(&qubits[i],&qubits[high]);
    return i;
}

void sort(BDDVAR* qubits, int low, int high)
{
    if (low < high) {
        int pi = partition(qubits, low, high);
        sort(qubits,low,pi-1);
        sort(qubits,pi+1,high);
    }
}

void handle_gate(C_struct* c_s, BDDVAR n_qubits, BDDVAR* qubits, Gate gate_s)
{
    // Initialise variables
    Gate gate_ctrl_s;
    BDDVAR j = 0, target = qubits[n_qubits];
    BDDVAR *control;
    // Increment the depth since were adding a column
    c_s->depth += 1;
    // Set the controlsize
    gate_s.controlSize = n_qubits;
    // If n_qubits is 0, it is not controlled
    if (n_qubits != 0) {
        // Create a list of control indices and store in <gate_s>
        control = malloc((n_qubits) * sizeof(BDDVAR));
        for (BDDVAR i = 0; i < n_qubits; i++) control[i] = qubits[i];
        gate_s.control = control;
        // sort qubits to place controls and in-betweens
        sort(qubits, 0, n_qubits);
        // Store all control gates in <c_s>
        for (BDDVAR i = qubits[0]; i <= qubits[n_qubits]; i++) {
            if (i == qubits[j]) {
                // Store target in control array of control gate (needed for optimisation)
                gate_ctrl_s = gate_ctrl;
                gate_ctrl_s.control = malloc(sizeof(BDDVAR));
                gate_ctrl_s.control[0] = target;
                gate_ctrl_s.controlSize = 1;
                // Store control gate in <c_s>
                c_s->circuit[i][c_s->depth] = gate_ctrl_s;
                // Increment j to check for next index
                j++;
            }
            // Once you find a control qubit, place vertical bars between controls and target
            else
                c_s->circuit[i][c_s->depth] = gate_ctrl_c;
        }
    }
    // Store gate on target qubit in <c_s>
    c_s->circuit[target][c_s->depth] = gate_s;
}

BDDVAR get_next_gate(C_struct* c_s, BDDVAR q, BDDVAR depth, bool successive)
{
    // Walk over the wire until a gate is found that is not a barrier or identity gate
    do {
        // If you reach the front or end of the circuit, a gate cannot be found
        if (depth == 0 || depth >= c_s->depth) return c_s->depth;
        depth += 2*successive-1;
    }
    while (c_s->circuit[q][depth].id == gate_barrier.id || c_s->circuit[q][depth].id == gate_I.id);
    return depth;
}

void optimize_c_struct(C_struct* c_s)
{
    // Initialise variables
    BDDVAR depth2;
    Gate gate;
    // Loop over all places in the circuit
    for (BDDVAR q = 0; q < c_s->qubits; q++) {
        for (BDDVAR depth1 = 0; depth1 < c_s->depth; depth1++) {
            // No need to optimize identity gates or control gates (controls are done along with its target)
            gate = c_s->circuit[q][depth1];
            if (gate.id != gate_I.id && gate.id != gate_ctrl.id && gate.id != gate_ctrl_c.id) {
                // Get the depth of the successive gate compared to the gate at the first depth
                depth2 = get_next_gate(c_s, q, depth1, true);
                if (depth2 == c_s->depth) break;
                gate = c_s->circuit[q][depth2];
                // Check if these gates are the same, and if so, recursively optimize (since we dont know palindrome length)
                if (find_palindromes(c_s, q, depth1, depth2)) {
                    for (BDDVAR i = 0; i < gate.controlSize; i++)
                        optimize_c_struct_p(c_s, gate.control[i], depth1, depth2);
                    optimize_c_struct_p(c_s, q, depth1, depth2);
                }
            }
        }
    }
}

void optimize_c_struct_p(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    // Initialise variables
    Gate gate;
    // First remove the two negating gates found in the previous iteration
    remove_gate(c_s, q, depth1);
    remove_gate(c_s, q, depth2);
    // Find the first preceding gate on the qubit
    depth1 = get_next_gate(c_s, q, depth1, false);
    if (depth1 == c_s->depth) return;
    // Find the first successive gate on the qubit
    depth2 = get_next_gate(c_s, q, depth2, true);
    if (depth2 == c_s->depth) return;
    // If two gates have been found, check all qubits
    for (BDDVAR q = 0; q < c_s->qubits; q++) {
        gate = c_s->circuit[q][depth1];
        if (gate.id != gate_I.id && gate.id != gate_ctrl.id && gate.id != gate_ctrl_c.id) {
            // Check if these gates are the same, and if so, recursively optimize (since we dont know palindrome length)
            if(find_palindromes(c_s, q, depth1, depth2)) {
                for (BDDVAR i = 0; i < gate.controlSize; i++)
                    optimize_c_struct_p(c_s, gate.control[i], depth1, depth2);
                optimize_c_struct_p(c_s, q, depth1, depth2);
            }
        }
    }
}

bool find_palindromes(C_struct* c_s, BDDVAR q, BDDVAR depth1, BDDVAR depth2)
{
    // Initialise variables
    Gate gate1 = c_s->circuit[q][depth1];
    Gate gate2 = c_s->circuit[q][depth2];
    // If the two gates do not fit the negation criteria, return
    if (gate1.id != gate2.id || (gate1.rotation + gate2.rotation != 1 && gate1.rotation != -gate2.rotation))
        return false;
    // Both gates must have the same controlsize, if not return false
    if(gate1.controlSize != gate2.controlSize)
        return false;
    // If the controlsize is 0, it is a single qubit gate and thus a negation
    if (gate1.controlSize == 0)
        return true;
    // Else also check all control qubits
    else {
        // Loop over all qubits in the control array
        for (BDDVAR i = 0; i < gate1.controlSize; i++) {
            // If the control qubits dont match, return false
            if (gate1.control[i] != gate2.control[i])
                return false;
            // Check all positions between these control gates
            for(BDDVAR j = depth1+1; j < depth2; j++) {
                // If they are not identity, return false (note that barriers should block optimises)
                if(c_s->circuit[gate1.control[i]][j].id != gate_I.id)
                    return false;
            }
        }
    }
    // If all tests above succeed and they havent returned false, the gates must be negations
    return true;
}

void remove_gate(C_struct* c_s, BDDVAR q, BDDVAR depth)
{
    // Remove all controls by placing identity gates
    if (c_s->circuit[q][depth].controlSize != 0) {
        for (BDDVAR i = c_s->circuit[q][depth].control[0]; i < q; i++) {
            c_s->circuit[i][depth] = gate_I;
        }
    }
    // Remove the gate by placing identity gates
    c_s->circuit[q][depth] = gate_I;
}

void reduce_c_struct(C_struct* c_s)
{
    // Initialise variables
    BDDVAR k = 0;
    // Reduce depth by moving gates to the left if possible
    for (BDDVAR j = 1; j <= c_s->depth; j++) {
        for (BDDVAR i = 0; i < c_s->qubits; i++) {
            if ((c_s->circuit[i][j].id != gate_I.id && c_s->circuit[i][j].id != gate_ctrl.id && 
            c_s->circuit[i][j].id != gate_ctrl_c.id && c_s->circuit[i][j].id != gate_barrier.id) || 
            (c_s->circuit[i][j].id == gate_barrier.id && i == c_s->qubits-1)) {
                if (c_s->circuit[i][j].classical_expect == -1)
                    reduce_gate(c_s, i, j);
                else
                    reduce_classical_gate(c_s, i, j);
            }
        }
    }
    // Get new depth
    for (BDDVAR i = 0; i < c_s->qubits; i++) {
        for (BDDVAR j = c_s->depth; j-- > 0;) {
            if (c_s->circuit[i][j].id != gate_I.id) {
                k = j < k ? k : j;
                break;
            }
        }
    }
    // Set new depth
    c_s->depth = k+1;
}

void reduce_classical_gate(C_struct* c_s, BDDVAR target, BDDVAR depth)
{
    // Initialise variables
    BDDVAR curr, reduce = depth, min = target, max = target;
    BDDVAR *control = c_s->circuit[target][depth].control;
    // Get the new reducable depth
    for (BDDVAR i = 0; i < c_s->qubits; i++) {
        curr = depth;
        // While possible, reduce depth by one
        while (curr > 0 && (c_s->circuit[i][curr-1].id == gate_I.id))
            curr--;
        if (curr > reduce) reduce = curr;
    }
    reduce--;
    // If the controls and gate are not reducable, return
    if (reduce - depth <= 0) return;
    // If the gate has controls, also process those
    if (c_s->circuit[target][depth].controlSize != 0) {
        // Check if target is below controls
        for (BDDVAR k = 0; k < c_s->circuit[target][depth].controlSize; k++) {
            min = (control[k] < min) ? control[k] : min;
            max = (control[k] > max) ? control[k] : max;
        }
        // Reduce the depth of all gates between first control and target
        for (BDDVAR i = min; i < max+1; i++) {
            if (c_s->circuit[i][depth].id != gate_I.id && i != target) {
                c_s->circuit[i][reduce] = c_s->circuit[i][depth];
                c_s->circuit[i][depth] = gate_I;
            }
        }
    }
    // Reduce the depth of the gate
    c_s->circuit[target][reduce] = c_s->circuit[target][depth];
    c_s->circuit[target][depth] = gate_I;
}

void reduce_gate(C_struct* c_s, BDDVAR target, BDDVAR depth)
{
    // Initialise variables
    BDDVAR curr, min = target, max = target;
    BDDVAR *control = c_s->circuit[target][depth].control;
    // Get the new reducable depth
    BDDVAR reduce = get_reduce_depth(c_s, target, depth);
    // If the gate has controls, also process those
    if (c_s->circuit[target][depth].controlSize != 0) {
        // Check if target is below controls
        for (BDDVAR k = 0; k < c_s->circuit[target][depth].controlSize; k++) {
            min = (control[k] < min) ? control[k] : min;
            max = (control[k] > max) ? control[k] : max;
        }
        // Check everything between first control and target
        for (BDDVAR i = min; i < max+1; i++) {
            // If the new reducable depth is less than the current, set current to new
            curr = get_reduce_depth(c_s, i, depth);
            if (curr > reduce) reduce = curr;
        }
        // If the controls and gate are not reducable, return
        if (reduce - depth <= 0) return;
        // Reduce the depth of all gates between first control and target
        for (BDDVAR i = min; i < max+1; i++) {
            if (c_s->circuit[i][depth].id != gate_I.id && i != target) {
                c_s->circuit[i][reduce] = c_s->circuit[i][depth];
                c_s->circuit[i][depth] = gate_I;
            }
        }
    }
    // If the gate is not reducable, return
    if (reduce - depth <= 0) return;
    // Reduce the depth of the gate
    c_s->circuit[target][reduce] = c_s->circuit[target][depth];
    c_s->circuit[target][depth] = gate_I;
}

BDDVAR get_reduce_depth(C_struct* c_s, BDDVAR target, BDDVAR depth)
{
    // While possible, reduce depth by one
    while (depth > 0 && (c_s->circuit[target][depth-1].id == gate_I.id))
        depth--;
    return depth;
}

void print_c_struct(C_struct c_s, bool show_rotation)
{
    // Initialise variables
    bool has_rotation = false;
    bool negative_rotation = false;
    // Loop over all positions, going row by row
    for (BDDVAR i = 0; i < c_s.qubits; i++) {
        for (BDDVAR j = 0; j < c_s.depth; j++) {
            // Print the gate
            if (c_s.circuit[i][j].classical_expect != -1)
                printf("c%s",c_s.circuit[i][j].gateSymbol);
            else
                printf("-%s",c_s.circuit[i][j].gateSymbol);
            // If rotation, check if a rotation needs to be printed in this column
            if (show_rotation) {
                // Reset rotation variables
                has_rotation = false;
                negative_rotation = false;
                // Check if any qubit in the column has a rotation and set variables accordingly
                for (BDDVAR k = 0; k < c_s.qubits; k++) {
                    if(c_s.circuit[k][j].id == gate_Rx.id || c_s.circuit[k][j].id == gate_Ry.id || c_s.circuit[k][j].id == gate_Rz.id)
                        has_rotation = true;
                    if(c_s.circuit[k][j].rotation < 0)
                        negative_rotation = true;
                }
                // If a rotation has been found in the column, print accordingly
                if (has_rotation) {
                    // If the current gate has a rotation, print its rotation
                    if(c_s.circuit[i][j].id == gate_Rx.id || c_s.circuit[i][j].id == gate_Ry.id || c_s.circuit[i][j].id == gate_Rz.id) {
                        printf("(%.4lf)",roundf(c_s.circuit[i][j].rotation*10000)/10000);
                        if(negative_rotation && c_s.circuit[i][j].rotation >= 0)
                            printf("-");
                    }
                    // If the current gate does not have a rotation, print extra wire space for alignment
                    else {
                        if (negative_rotation)
                            printf("---------");
                        else
                            printf("--------");
                    }
                }
            }
        }
        // After a wire has been printed, go to newline
        printf("-\n");
    }
}
