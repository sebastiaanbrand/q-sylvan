#include "grover_cnf.h"
#include "sylvan_qdd_complex.h"

TASK_IMPL_5(QDD, qdd_grover_cnf_iteration, QDD, qdd, BDDVAR, n, BDDVAR, k, BDDVAR, clauses, int*, oracle)
{
    char *filename = malloc(strlen("../plotting/Grover_.txt")+10);
    sprintf(filename, "../plotting/Grover%d_%d.txt", n, clauses);
    FILE *f = fopen(filename, "a");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    // Compute the results of each clause
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[clause*k+l] < 0) {
                qdd = qdd_gate(qdd, GATEID_X, abs(oracle[clause*k+l])-1);
                fprintf(f, "%d: X\n",abs(oracle[clause*k+l])-1);
            }
        }
        // "oracle" call (Save results of clause in ancilla)
        qdd = qdd_cgate3(qdd,GATEID_X,abs(oracle[clause*k])-1,abs(oracle[clause*k+1])-1,abs(oracle[clause*k+2])-1,n+1+clause);
        fprintf(f, "%d/%d/%d->%d: X\n",abs(oracle[clause*k])-1,abs(oracle[clause*k+1])-1,abs(oracle[clause*k+2])-1,n+1+clause);
        // Uncompute X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[clause*k+l] < 0) {
                qdd = qdd_gate(qdd, GATEID_X, abs(oracle[clause*k+l])-1);
                fprintf(f, "%d: X\n",abs(oracle[clause*k+l])-1);
            }
        }
        fprintf(f, "barrier\n");
    }
    // Oracle flip (since the MCX is in the wrong direction, we use a Z and 2*H around the target)
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        qdd = qdd_gate(qdd, GATEID_X, n+1+clause);
        fprintf(f, "%d: X\n",n+1+clause);
    }
    qdd = qdd_gate(qdd, GATEID_H, n);
    fprintf(f, "%d: H\n",n);
    qdd = qdd_cgate_range(qdd,GATEID_Z,n,n+clauses-1,n+clauses);
    for(BDDVAR i = n; i < n+clauses-1; i++) { fprintf(f, "%d/",i); }
    fprintf(f, "%d->%d: Z\n",n+clauses-1,n+clauses);
    qdd = qdd_gate(qdd, GATEID_H, n);
    fprintf(f, "%d: H\n",n);
    for (BDDVAR clause = 0; clause < clauses; clause++) {
        qdd = qdd_gate(qdd, GATEID_X, n+1+clause);
        fprintf(f, "%d: X\n",n+1+clause);
    }
    fprintf(f, "barrier\n");
    // Uncompute all clauses
    for (BDDVAR clause = clauses; clause != 0; clause--) {
        // Satisfy clause rules by applying X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[(clause-1)*k+l] < 0) {
                qdd = qdd_gate(qdd, GATEID_X, abs(oracle[(clause-1)*k+l])-1);
                fprintf(f, "%d: X\n",abs(oracle[(clause-1)*k+l])-1);
            }
        }
        qdd = qdd_cgate3(qdd,GATEID_X,abs(oracle[(clause-1)*k])-1,abs(oracle[(clause-1)*k+1])-1,abs(oracle[(clause-1)*k+2])-1,n+1+(clause-1));
        fprintf(f, "%d/%d/%d->%d: X\n",abs(oracle[(clause-1)*k])-1,abs(oracle[(clause-1)*k+1])-1,abs(oracle[(clause-1)*k+2])-1,n+1+(clause-1));
        // Uncompute X gates
        for (BDDVAR l = 0; l < k; l++) {
            if(oracle[(clause-1)*k+l] < 0) {
                qdd = qdd_gate(qdd, GATEID_X, abs(oracle[(clause-1)*k+l])-1);
                fprintf(f, "%d: X\n",abs(oracle[(clause-1)*k+l])-1);
            }
        }
        fprintf(f, "barrier\n");
    }

    // H on all qubits (also ancilla to make CZ a CX)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qdd = qdd_gate(qdd, GATEID_H, qubit);
        fprintf(f, "%d: H\n",qubit);
    }

    // X on all qubits (except ancilla)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qdd = qdd_gate(qdd, GATEID_X, qubit);
        fprintf(f, "%d: X\n",qubit);
    }

    // Controlled-NOT over all qubits (on ancilla)
    qdd = qdd_cgate_range(qdd, GATEID_X, 0, n-1, n);
    for(BDDVAR i = 0; i < n-1; i++) { fprintf(f, "%d/",i); }
    fprintf(f, "%d->%d: X\n",n-1,n);

    // X on all qubits (except ancilla)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qdd = qdd_gate(qdd, GATEID_X, qubit);
        fprintf(f, "%d: X\n",qubit);
    }

    // H on all qubits (also ancilla to maqubite CZ a CX)
    for (BDDVAR qubit = 0; qubit < n; qubit++) {
        qdd = qdd_gate(qdd, GATEID_H, qubit);
        fprintf(f, "%d: H\n",qubit);
    }
    fprintf(f, "barrier\n");
    fclose(f);

    return qdd;
}

QDD
qdd_grover_cnf(BDDVAR n, int* oracle, BDDVAR k, BDDVAR clauses, BDDVAR n_answers)
{
    char *filename = malloc(strlen("../plotting/Grover_.txt")+10);
    sprintf(filename, "../plotting/Grover%d_%d.txt", n, clauses);
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }
    LACE_ME;
    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) / n_answers ) );
    // start with all zero state + ancilla: |000...0>|1>
    bool *init = malloc( sizeof(bool)*(n+1+clauses) );
    for (BDDVAR qubit = 0; qubit < n+1+clauses; qubit++) init[qubit] = 0;
    // init[n] = 1; 
    fprintf(f, "%d qubits\n", n+1+clauses);
    QDD qdd = qdd_create_basis_state(n+1+clauses, init);
    qdd = qdd_gate(qdd, GATEID_X, n);
        fprintf(f, "%d: X\n", n);

    // H on all qubits
    for (BDDVAR qubit = 0; qubit < n+1; qubit++) {
        qdd = qdd_gate(qdd, GATEID_H, qubit);
        fprintf(f, "%d: H\n", qubit);
    }
    fprintf(f, "barrier\n");
    fclose(f);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_cnf_iteration(qdd, n, k, clauses, oracle);
    }

    return qdd;
}
