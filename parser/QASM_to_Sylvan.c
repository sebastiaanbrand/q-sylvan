#include "QASM_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"


// void remove_multiline_comments(char *str)
// {
//     char *c = str;
//     c=strchr(c,'/');
//     if (c != NULL) {
//         if (str[c-str+1] == '*') {
//         }
//     }
// }

TASK_IMPL_4(int, final_measuring, QDD, qdd, int16_t*, measurements, BDDVAR*, nvars, BDDVAR, shots)
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
    for (BDDVAR i = 0; i < *nvars; i++) { if (measurements[i] != -1) sum++; }
    int *probs = malloc(pow(2, sum) * sizeof(int));
    for (int k = 0; k < (1 << (*nvars)); k++) {
        bool *x = int_to_bitarray(k, *nvars, false);
        j = sum-1;
        index = 0;
        for (BDDVAR i = *nvars; i > 0; i--)
            if (measurements[i-1] != -1) { index += x[i-1] * pow(2, j--); }
        probs[index] += hits[k];
    }

    // Print reformatted circuit results
    for(int k = 0; k < pow(2, sum); k++) {
        bool *x = int_to_bitarray(k, *nvars, false);
        j = sum-1;
        for (BDDVAR i = *nvars; i > 0 ; i--) {
            if (measurements[i-1] != -1) {
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

TASK_IMPL_4(QDD, handle_intermediate_measure, QDD, qdd, int16_t*, measurements, BDDVAR, nvars, bool*, creg)
{
    int *m = malloc(sizeof(int));
    double *p = malloc(sizeof(double));
    for (BDDVAR i = 0; i < nvars; i++) {
        if (measurements[i] != -1) {
            qdd = qdd_measure_qubit(qdd, i, nvars, m, p);
            creg[measurements[i]] = *m;
            measurements[i] = -1;
        }
    }
    return qdd;
}

TASK_IMPL_5(QDD, handle_line, QDD, qdd, char*, line, int16_t*, measurements, BDDVAR*, nvars, bool*, creg)
{
    uint32_t *gateid = malloc(sizeof(uint32_t));
    *gateid = 0;
    bool isgate = false;
    BDDVAR n_qubits = 0;
    char *temp;
    char *tokens[2];
    // Create all-zero state with "temp" qubits
    if (strstr(line, "qreg") != NULL) {
        get_qubits(line, 0, nvars);
        creg = malloc(sizeof(bool)*(*nvars));
        for(BDDVAR i = 0; i < *nvars; i++)
            creg[i] = 0;
        qdd = qdd_create_all_zero_state(*nvars);
    }
    // Set measure marker on qubit "temp"
    else if (strstr(line, "measure") != NULL) {
        BDDVAR *measure = malloc(2*sizeof(BDDVAR));
        get_qubits(line, 1, measure);
        measurements[measure[0]] = (int16_t)measure[1];
    }
    // Handle if statement
    else if (strstr(line, "if") != NULL) {
        qdd = handle_intermediate_measure(qdd, measurements, *nvars, creg);
        BDDVAR *sum = malloc(sizeof(BDDVAR));
        *sum = 0;
        temp = strtok(line, "=");
        line = strtok(NULL, "");
        if(strstr(temp, "[") != NULL) {
            get_qubits(temp, 0, sum);
            *sum = creg[*sum];
        }
        else {
            for (BDDVAR i = 0; i < *nvars; i++)
                *sum += creg[i]*pow(2,i);
        }
        temp = strtok(line, "=");
        temp = strtok(temp, ")");
        if (atoi(temp) == (int)*sum) {
            line = strtok(NULL, "");
            goto handleGate;
        }
    }
    // // Handle for loop
    // else if (strstr(line, "for") == 0) {
    //     temp = strtok(line, "[");
    //     temp = strtok(NULL, "]");
    //     measurements[atoi(temp)] = true;
    // }
    else {
        handleGate:
        // tokenize string
        tokens[0] = strtok(line, " ");
        tokens[1] = strtok(NULL, "");
        isgate = get_gateid(tokens[0], gateid);
        if (!isgate)
            return qdd;
        while(tokens[0][n_qubits] == 'c' && strcmp(tokens[0], "creg") != 0)
            n_qubits++;
        BDDVAR *qubits = malloc((n_qubits+1) * sizeof(BDDVAR));
        get_qubits(tokens[1], n_qubits, qubits);

        qdd = handle_intermediate_measure(qdd, measurements, *nvars, creg);

        if (n_qubits == 0)
            qdd = qdd_gate(qdd, *gateid, qubits[0]);
        else if (n_qubits == 1)
            qdd = qdd_cgate(qdd, *gateid, qubits[0], qubits[1]);
        else if (n_qubits == 2)
            qdd = qdd_cgate2(qdd, *gateid, qubits[0], qubits[1], qubits[2]);
        else if (n_qubits == 3) {
            qdd = qdd_cgate3(qdd, *gateid, qubits[0], qubits[1], qubits[2], qubits[3]);
        }
        else exit(-1);
    }
    return qdd;
}

TASK_IMPL_3(bool, get_qubits, char*, token, BDDVAR, n_qubits, BDDVAR*, qubits)
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

TASK_IMPL_2(bool, get_gateid, char*, gate, uint32_t*, gate_id)
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

void read_QASM(char *filename, BDDVAR shots)
{
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        perror("Error while opening the file.\n");
        exit(-1);
    }

    LACE_ME;

    char *line = NULL, *c;
    BDDVAR *nvars = malloc(sizeof(BDDVAR));
    *nvars = 0;
    bool *creg = malloc(1024 * sizeof(bool));

    // Since we dont know yet how many qubits the circuit will contain,
    // we initialize measurements with size 1024 (since bool takes up minimal space)
    int16_t *measurements = malloc(1024 * sizeof(int16_t));
    for (BDDVAR i = 0; i < 1024; i++)
    {
        measurements[i] = -1;
    }
    size_t len = 0;
    ssize_t read;
    QDD qdd = qdd_create_all_zero_state(0);
    while ((read = getline(&line, &len, f)) != -1) {
        // skip if comment
        if(line[0] == '/' && line[1] == '/')
            continue;
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t'))
            line++;
        // remove empty lines, trailing information after ';'
        c = strchr(line, ';');
        if (c != NULL)
        {
            line[c - line] = '\0';
            qdd = handle_line(qdd, line, measurements, nvars, creg);
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
    BDDVAR seed = 100;
    int opt;
    while((opt = getopt(argc, argv, "f:s:r:")) != -1) {
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

    read_QASM(filename, shots);

    sylvan_quit();
    lace_exit();

    return 0;
}