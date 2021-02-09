#include "QASM_to_circuit.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

TASK_IMPL_3(bool, get_qubits_circuit, char*, token, BDDVAR, n_qubits, BDDVAR*, qubits)
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

TASK_IMPL_2(bool, get_gateid_circuit, char*, gate, uint32_t*, gate_id)
{
    while(gate[0] == 'c' && strcmp(gate, "creg") != 0) {
        gate++;
    }
    // Handle gates
    if (strcmp(gate, "i") == 0)
        *gate_id = GATEID_I;
    else if (strcmp(gate, "h") == 0)
        *gate_id = GATEID_H;
    else if (strcmp(gate, "x") == 0)
        *gate_id = GATEID_X;
    else if (strcmp(gate, "y") == 0)
        *gate_id = GATEID_Y;
    else if (strcmp(gate, "z") == 0)
        *gate_id = GATEID_Z;
    else if (strcmp(gate, "s") == 0)
        *gate_id = GATEID_S;
    else if (strcmp(gate, "sdg") == 0)
        *gate_id = GATEID_Rk_dag(-2);
    else if (strcmp(gate, "t") == 0)
        *gate_id = GATEID_T;
    else if (strcmp(gate, "tdg") == 0)
        *gate_id = GATEID_Tdag;
    else {
        gate = strtok(gate, "(");
        if (strcmp(gate, "rx") == 0)
            *gate_id = GATEID_Rx(atof(strtok(NULL, ")")));
        else if (strcmp(gate, "ry") == 0)
            *gate_id = GATEID_Ry(atof(strtok(NULL, ")")));
        else if (strcmp(gate, "rz") == 0)
            *gate_id = GATEID_Rz(atof(strtok(NULL, ")")));
        else {
            return false;
        }
    }
    return true;
}

TASK_IMPL_5(bool, get_parallel_depth, char*, targets, int32_t**, circuit, BDDVAR, n_qubits, BDDVAR*, curr_depth, uint32_t*, gateid)
{
    bool parallelize = true;
    BDDVAR *qubits = malloc((n_qubits+1) * sizeof(BDDVAR));
    get_qubits_circuit(targets, n_qubits, qubits);
    for (BDDVAR i = 0; i <= n_qubits; i++) {
        if (circuit[qubits[i]][*curr_depth] != 0) {
            parallelize = false;
            *curr_depth += 1;
            break;
        }
    }
    BDDVAR depth = *curr_depth;
    while (parallelize && depth > 0) {
        depth--;
        for (BDDVAR i = 0; i <= n_qubits; i++) {
            if (circuit[qubits[i]][depth] != 0) {
                parallelize = false;
                depth++;
                break;
            }
        }
    }
    for (BDDVAR i = 0; i < n_qubits; i++)
        circuit[qubits[i]][depth] = -qubits[n_qubits];
    circuit[qubits[n_qubits]][depth] = *gateid;
    return true;
}

TASK_IMPL_4(bool, handle_line_circuit, char*, line, int32_t**, circuit, BDDVAR*, nvars, BDDVAR*, curr_depth)
{
    uint32_t *gateid = malloc(sizeof(uint32_t));
    *gateid = -*nvars;
    bool isgate = false;
    BDDVAR n_qubits = 0;
    char *tokens[2];
    // tokenize string
    tokens[0] = strtok(line, " ");
    tokens[1] = strtok(NULL, "");
    // Create all-zero state with "temp" qubits
    if (strstr(tokens[0], "qreg") != NULL) {
        get_qubits_circuit(tokens[1], 0, nvars);
    }
    // Set measure marker on qubit "temp"
    else if (strstr(line, "measure") != NULL) {
        get_parallel_depth(tokens[1], circuit, n_qubits, curr_depth, gateid);
    }
    else {
        isgate = get_gateid_circuit(tokens[0], gateid);
        if (!isgate)
            return false;
        while(tokens[0][n_qubits] == 'c' && strcmp(tokens[0], "creg") != 0)
            n_qubits++;
        get_parallel_depth(tokens[1], circuit, n_qubits, curr_depth, gateid);
    }
    return true;
}

void make_circuit(char *filename, bool optimize)
{
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        perror("Error while opening the file.\n");
        return;
    }

    LACE_ME;

    char *line = NULL, *c;
    BDDVAR max_qubits = 128;
    BDDVAR max_wire = 1024;
    int32_t **circuit = malloc(max_qubits*sizeof(*circuit));
    for (BDDVAR i = 0; i < max_qubits; ++i) {
        circuit[i] = malloc(max_wire * sizeof *circuit[i]);
        for (BDDVAR j = 0; j < max_wire; j++)
            circuit[i][j] = 0; // GATEID_I
    }
    BDDVAR *nvars = malloc(sizeof(BDDVAR));
    BDDVAR *curr_depth = malloc(sizeof(BDDVAR));
    *nvars = 0;
    *curr_depth = 0;
    size_t len = 0;
    ssize_t read;
    while ((read = getline(&line, &len, f)) != -1) {
        if(*nvars >= max_qubits) {
            printf("Maximum number of qubits is 128. Currently exceeding.");
            return;
        }
        if(*curr_depth >= max_wire) {
            printf("Maximum wire depth is 1024. Currently exceeding.");
            return;
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
    print_circuit(circuit, nvars, curr_depth);
}

BDDVAR optimize_circuit(int32_t** circuit, BDDVAR nvars, BDDVAR depth)
{
    for (BDDVAR j = 0; j < depth; j++) {
        for (BDDVAR i = 0; i < nvars; i++) {
            if (circuit[i][j] != 0)
                find_palindromes(circuit, i, j, nvars, depth);
        }
    }
    return reduce_circuit(circuit, nvars, depth);
}

void find_palindromes(int32_t** circuit, BDDVAR first_q, BDDVAR curr_depth, BDDVAR nvars, BDDVAR depth)
{
    if (circuit[first_q][curr_depth] == 0) return;
    bool *qubits = malloc(nvars*sizeof(bool));
    BDDVAR end_depth, j;
    bool found_palindrome = false;
    end_depth = curr_depth;
    loop:
    for (BDDVAR i = 0; i < nvars; i++) qubits[i] = false;
    qubits[first_q] = true;
    j = end_depth+1;
    while(j != depth && !found_palindrome) {
        j = end_depth+1;
        while (j != depth && circuit[first_q][j] != circuit[first_q][curr_depth]) j++;
        if (j == depth) return;
        end_depth = j;
        for (BDDVAR i = curr_depth; i < j; i++,j--) {
            for (BDDVAR q = 0; q < nvars; q++) {
                if (circuit[q][i] < 0 && circuit[q][i] != -(int32_t)nvars)
                    qubits[-circuit[q][i]] = true;
                if (circuit[q][j] < 0 && circuit[q][j] != -(int32_t)nvars)
                    qubits[-circuit[q][j]] = true;
                if ((circuit[q][i] < 0 && qubits[-circuit[q][i]] && circuit[q][i] != -(int32_t)nvars) ||
                        (circuit[q][j] < 0 && qubits[-circuit[q][j]] && circuit[q][j] != -(int32_t)nvars))
                    qubits[q] = true;
                if (qubits[q] && circuit[q][i] != circuit[q][j])
                    goto loop;
            }
        }
        found_palindrome = true;
    }
    for (BDDVAR i = 0; i < nvars; i++) qubits[i] = false;
    qubits[first_q] = true;
    j = end_depth;
    for (BDDVAR i = curr_depth; i < j; i++,j--) {
        for (BDDVAR q = 0; q < nvars; q++) {
            if (circuit[q][i] < 0 && circuit[q][i] != -(int32_t)nvars)
                qubits[-circuit[q][i]] = true;
            if (circuit[q][j] < 0 && circuit[q][j] != -(int32_t)nvars)
                qubits[-circuit[q][j]] = true;
            if ((circuit[q][i] < 0 && qubits[-circuit[q][i]] && circuit[q][i] != -(int32_t)nvars) ||
                    (circuit[q][j] < 0 && qubits[-circuit[q][j]] && circuit[q][j] != -(int32_t)nvars))
                qubits[q] = true;
            if (qubits[q]) {
                circuit[q][i] = 0;
                circuit[q][j] = 0;
            }
        }
    }
}

BDDVAR reduce_circuit(int32_t** circuit, BDDVAR nvars, BDDVAR depth)
{
    int32_t k = 0;
    int32_t is_controlled;
    for (BDDVAR j = 1; j <= depth; j++) {
        is_controlled = -1;
        for (BDDVAR i = 0; i < nvars; i++) {
            if (circuit[i][j] < 0 && circuit[i][j] != -(int32_t)nvars) {
                reduce_controlled_gate(circuit, j, -circuit[i][j], nvars);
                is_controlled = -circuit[i][j];
            }
            else if (circuit[i][j] != 0) {
                k = j;
                while (circuit[i][k-1] == 0 && k > 0) k--;
                if ((int32_t)i != is_controlled) {
                    circuit[i][k] = circuit[i][j];
                    circuit[i][j] = 0;
                }
            }
        }
    }
    k = 0;
    for (BDDVAR i = 0; i < nvars; i++) {
        for (BDDVAR j = depth; j-- > 0;) {
            if (circuit[i][j] != 0) {
                k = (int32_t)j < k ? k : (int32_t)j;
                break;
            }
        }
    }
    return k;
}

void reduce_controlled_gate(int32_t** circuit, BDDVAR depth, BDDVAR target, BDDVAR nvars)
{
    BDDVAR j = depth;
    for (; j-- > 0;) {
        for (BDDVAR i = 0; i < nvars; i++) {
            if((circuit[i][depth] < 0 || i == target) && circuit[i][j] != 0)
                goto reduce;
        }
    }
    reduce:
    j++;
    if (j == depth)
        return;
    else {
        for (BDDVAR i = 0; i < nvars; i++) {
            if(circuit[i][depth] < 0 || i == target) {
                circuit[i][j] = circuit[i][depth];
                circuit[i][depth] = 0;
            }
        }
    }
}

void print_circuit(int32_t** circuit, BDDVAR* nvars, BDDVAR* curr_depth)
{
    for (BDDVAR i = 0; i < *nvars; i++) {
        for (BDDVAR j = 0; j < *curr_depth+1; j++) {
            if (circuit[i][j] >= 0 && circuit[i][j] < 9)
                printf("_");
            if ((circuit[i][j] >= 0 && circuit[i][j] < 99) || (circuit[i][j] > -9 && circuit[i][j] < 0))
                printf("_");
            printf("_%d",circuit[i][j]);
            if ((circuit[i][j] >= 0 && circuit[i][j] < 999) || (circuit[i][j] > -99 && circuit[i][j] < 0))
                printf("_");
        }
        printf("_\n");
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

    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<16, -1, true);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    make_circuit(filename, optimize);

    sylvan_quit();
    lace_exit();

    return 0;
}