#include "grover.h"
#include "sylvan_qdd_complex.h"


TASK_IMPL_3(QDD, qdd_grover_iteration, QDD, qdd, BDDVAR, n, bool*, flag)
{
    // "oracle" call  (apply -1 flag to desired amplitude)
    qdd = qdd_all_control_phase(qdd, n, flag);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Phase on all amplitudes except |000...0>
    bool x[n]; 
    for(BDDVAR k = 0; k < n; k++) x[k] = 0;
    qdd = qdd_all_control_phase(qdd, n, x);
    qdd = qdd_scalar_mult(qdd, comp_minus_one());

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    return qdd;
}

QDD
qdd_grover(BDDVAR n, bool* flag)
{   
    LACE_ME;

    // not entirely sure about this, book says R <= ceil(pi/4 * sqrt(N))
    uint32_t R = floor( 3.14159265359/4.0 * sqrt( pow(2,n) ) );

    // start with all zero state |000...0>
    QDD qdd = qdd_create_all_zero_state(n);

    // H on all qubits
    for (BDDVAR k = 0; k < n; k++) qdd = qdd_gate(qdd, GATEID_H, k);

    // Grover iterations
    for (uint32_t i = 1; i <= R; i++) {
        qdd = qdd_grover_iteration(qdd, n, flag);
    }

    return qdd;
}

bool *
qdd_grover_random_flag(BDDVAR nqubits)
{
    bool *flag = malloc( (sizeof(bool)*nqubits) );
    for (BDDVAR k = 0; k < nqubits; k++) flag[k] = (rand() % 2);
    return flag;
}
