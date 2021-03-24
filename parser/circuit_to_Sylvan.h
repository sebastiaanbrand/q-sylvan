#include <stdbool.h>

#include "QASM_to_circuit.h"
#include "sylvan.h"

#define get_gate_id(gate) (CALL(get_gate_id,gate));
TASK_DECL_1(BDDVAR, get_gate_id, Gate);
#define apply_gate(qdd,gate,i) (CALL(apply_gate,qdd,gate,i));
TASK_DECL_3(QDD, apply_gate, QDD, Gate, BDDVAR);
#define handle_control(c_s,k,depth,gate_id) (CALL(handle_control,c_s,k,depth,gate_id));
TASK_DECL_4(QDD, handle_control, C_struct, BDDVAR, BDDVAR, BDDVAR);


bool measure(QDD qdd, bool* measurements, BDDVAR nvars, BDDVAR shots, bool show);
bool run_c_struct_matrix(C_struct c_s, BDDVAR shots, bool show);
bool run_c_struct(C_struct c_s, BDDVAR shots, bool show);
