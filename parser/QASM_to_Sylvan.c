#include "QASM_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// void remove_multiline_comments(char *str)
// {
//     char *c = str;
//     c=strchr(c,'/');
//     if (c != NULL) {
//         if (str[c-str+1] == '*') {
//         }
//     }
// }

TASK_IMPL_4(int, final_measuring, QDD, qdd, bool*, measurements, BDDVAR*, nvars, BDDVAR, shots)
{
    // Run the circuit x times, where x == shots
    bool *ms = malloc(*nvars * sizeof(bool));
    double *p = malloc(*nvars * sizeof(double));
    int *hits = malloc(pow(2,*nvars) * sizeof(int));
    for(int i = 0; i < pow(2,*nvars); i++) { hits[i] = 0; }
    for(BDDVAR i = 0; i < shots; i++) { 
        qdd_measure_all(qdd, *nvars, ms, p);
        hits[bitarray_to_int(ms, *nvars, false)]++;
    }

    // Reformat circuit results based on what qubits were measured
    BDDVAR index, j, sum = 0;
    for (BDDVAR i = 0; i < *nvars; i++) { if (measurements[i]) sum++; }
    int *probs = malloc(pow(2, sum) * sizeof(int));
    for (int k = 0; k < (1 << (*nvars)); k++) {
        bool *x = int_to_bitarray(k, *nvars, false);
        j = sum-1;
        index = 0;
        for (BDDVAR i = *nvars; i > 0; i--)
            if (measurements[i-1]) { index += x[i-1] * pow(2, j--); }
        probs[index] += hits[k];
    }

    // Print reformatted circuit results
    for(int k = 0; k < pow(2, sum); k++) {
        bool *x = int_to_bitarray(k, *nvars, false);
        j = sum-1;
        for (BDDVAR i = *nvars; i > 0 ; i--) {
            if (measurements[i-1]) {
                printf("%d", x[j]);
                j--;
            }
            else
                printf("_");
        }
        printf(": %d\n", probs[k]);
    }
    return 0;
}

TASK_IMPL_3(QDD, handle_intermediate_measure, QDD, qdd, bool*, measurements, BDDVAR, nvars)
{
    int *m = malloc(sizeof(int));
    double *p = malloc(sizeof(double));
    for (BDDVAR i = 0; i < nvars; i++)
    {
        if (measurements[i])
        {
            qdd = qdd_measure_qubit(qdd, i, nvars, m, p);
            measurements[i] = false;
        }
    }
    return qdd;
}

TASK_IMPL_3(QDD, handle_single_qubit_gate, QDD, qdd, char*, target, uint32_t, gate_id)
{
    target = strtok(target, "[");
    target = strtok(NULL, "]");
    return qdd_gate(qdd, gate_id, atoi(target));
}

TASK_IMPL_2(QDD, handle_controlled_qubit_gate, QDD, qdd, char**, tokens)
{
    BDDVAR i = 0;
    while(tokens[0][i] == 'c') {
        i++;
    }
    char **qubits = malloc((i+1) * sizeof(char));
    qubits[0] = strtok(tokens[1], "[");
    BDDVAR j = 0;
    for(; j < i; j++) {
        qubits[j] = strtok(NULL, "]");
        qubits[j+1] = strtok(NULL, "[");
    }
    qubits[j] = strtok(NULL, "]");
    if(i == 1)
        qdd = qdd_cgate(qdd, GATEID_X, atoi(qubits[0]), atoi(qubits[1]));
    else if(i == 2)
        qdd = qdd_cgate2(qdd, GATEID_X, atoi(qubits[0]), atoi(qubits[1]), atoi(qubits[2]));
    else if(i == 3)
        qdd = qdd_cgate3(qdd, GATEID_X, atoi(qubits[0]), atoi(qubits[1]), atoi(qubits[2]), atoi(qubits[3]));
    else exit(-1);
    return qdd;
}

TASK_IMPL_4(QDD, handle_tokens, QDD, qdd, char**, tokens, bool*, measurements, BDDVAR*, nvars)
{
    char *temp;
    // Create all-zero state with "temp" qubits
    if (strcmp(tokens[0], "qreg") == 0)
    {
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        *nvars = (BDDVAR)atoi(temp);
        qdd = qdd_create_all_zero_state(atoi(temp));
    }
    // Set measure marker on qubit "temp"
    else if (strcmp(tokens[0], "measure") == 0)
    {
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        measurements[atoi(temp)] = true;
    }
    // Handle gates
    if (strcmp(tokens[0], "h") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_H);
    }
    else if (strcmp(tokens[0], "x") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_X);
    }
    else if (strcmp(tokens[0], "y") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Y);
    }
    else if (strcmp(tokens[0], "z") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Z);
    }
    else if (strcmp(tokens[0], "s") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_S);
    }
    else if (strcmp(tokens[0], "sdg") == 0)
    {
        printf("TODO: sdg does not exist");
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        // qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID);
    }
    else if (strcmp(tokens[0], "t") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_T);
    }
    else if (strcmp(tokens[0], "tdg") == 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Tdag);
    }
    else if (strcmp(strtok(tokens[0], "("), "rx") == 0)
    {
        temp = strtok(NULL, ")");
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Rx(atof(temp)));
    }
    else if (strcmp(strtok(tokens[0], "("), "ry") == 0)
    {
        temp = strtok(NULL, ")");
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Ry(atof(temp)));
    }
    else if (strcmp(strtok(tokens[0], "("), "rz") == 0)
    {
        temp = strtok(NULL, ")");
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_single_qubit_gate(qdd, tokens[1], GATEID_Rz(atof(temp)));
    }
    else if (tokens[0][0] == 'c' && strcmp(tokens[0], "creg") != 0)
    {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars);
        qdd = handle_controlled_qubit_gate(qdd, tokens);
    }
    return qdd;
}

void read_QASM(char *filename, BDDVAR shots)
{
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL)
    {
        perror("Error while opening the file.\n");
        exit(-1);
    }

    LACE_ME;

    char *line = NULL, *c;
    BDDVAR *nvars = malloc(sizeof(BDDVAR));
    *nvars = 0;
    // Since we dont know yet how many qubits the circuit will contain,
    // we initialize measurements with size 1024 (since bool takes up minimal space)
    bool *measurements = malloc(1024 * sizeof(bool));
    for (BDDVAR i = 0; i < 1024; i++)
    {
        measurements[i] = false;
    }
    size_t len = 0;
    ssize_t read;
    char *tokens[2];
    QDD qdd = qdd_create_all_zero_state(0);
    while ((read = getline(&line, &len, f)) != -1)
    {
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t'))
            line++;
        // remove empty lines, trailing information after ';'
        c = strchr(line, ';');
        if (c != NULL)
        {
            line[c - line] = '\0';
            // tokenize string
            tokens[0] = strtok(line, " ");
            tokens[1] = strtok(NULL, "");
            qdd = handle_tokens(qdd, tokens, measurements, nvars);
        }
    }
    // test probabilities
    final_measuring(qdd, measurements, nvars, shots);
    fclose(f);
}

int main(int argc, char *argv[])
{
    // check for file
    char *filename = "";
    BDDVAR shots = 100;
    int opt;
    while((opt = getopt(argc, argv, "f:s:")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 's':
                shots = atoi(optarg);
                break;
        }
    }
    if(strcmp(filename, "") == 0)
    {
        printf("Give filename of qasm file.\n");
        exit(-1);
    }
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL << 25, 1LL << 25, 1LL << 16, 1LL << 16);
    sylvan_init_package();
    sylvan_init_qdd(1LL << 19, -1);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    read_QASM(filename, shots);

    sylvan_quit();
    lace_exit();

    return 0;
}