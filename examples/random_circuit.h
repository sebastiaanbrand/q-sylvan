#include <qsylvan.h>

void random_qubit(BDDVAR nqubits, BDDVAR *t);
void random_control_target(BDDVAR nqubits, BDDVAR *c, BDDVAR *t);
QDD qdd_run_random_circuit(BDDVAR nqubits, uint64_t ngates, double cgate_ratio, uint64_t rseed);
QDD qdd_run_random_single_qubit_gates(BDDVAR nqubits, uint64_t ngates, uint64_t rseed);