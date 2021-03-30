#include <stdbool.h>

#include "QASM_to_circuit.h"
#include "sylvan.h"

/**
 * Returns the Sylvan GATEID that corresponds to <gate>
 * 
 * @param gate Gate_struct of which to get the Sylvan GATEID
 * 
 * @return Sylvan GATEID which corresponds to <gate>
 */
#define get_gate_id(gate) (CALL(get_gate_id,gate));
TASK_DECL_1(BDDVAR, get_gate_id, Gate);

/**
 * Applies <gate> to <qdd> on the qubit with index <i>.
 * 
 * @param qdd the statevector on which to apply <gate>
 * @param gate the gate to be applied
 * @param i the index of the qubit in <qdd> on which to apply <gate>
 * 
 * @return the statevector qdd where <gate> has been applied on the qubit with index <i>
 */
#define apply_gate(qdd,gate,i) (CALL(apply_gate,qdd,gate,i));
TASK_DECL_3(QDD, apply_gate, QDD, Gate, BDDVAR);

#define handle_control(gate,k,n) (CALL(handle_control,gate,k,n));
TASK_DECL_3(QDD, handle_control, Gate, BDDVAR, BDDVAR);


void final_measure(QDD qdd, bool* measurements, BDDVAR nvars, BDDVAR shots, bool show);
QDD run_c_struct_matrix(C_struct c_s, BDDVAR shots, bool show);
QDD run_c_struct(C_struct c_s, BDDVAR shots, bool show);
