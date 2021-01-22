#include <stdbool.h>

#include "sylvan.h"

void read_QASM(char *filename, BDDVAR shots);

#define handle_intermediate_measure(qdd,measurements,nvars) (CALL(handle_intermediate_measure,qdd,measurements,nvars));
TASK_DECL_3(QDD, handle_intermediate_measure, QDD, bool*, BDDVAR);

#define handle_single_qubit_gate(qdd,target,gate_id) (CALL(handle_single_qubit_gate,qdd,target,gate_id));
TASK_DECL_3(QDD, handle_single_qubit_gate, QDD, char*, uint32_t);

#define get_gateid(tokens,gate_id) (CALL(get_gateid,tokens,gate_id));
TASK_DECL_2(uint32_t, get_gateid, char*, uint32_t*);

#define handle_tokens(qdd,tokens,measurements,nvars) (CALL(handle_tokens,qdd,tokens,measurements,nvars));
TASK_DECL_4(QDD, handle_tokens, QDD, char**,bool*, BDDVAR*);

#define final_measuring(qdd,measurements,nvars,shots) (CALL(final_measuring,qdd,measurements,nvars,shots));
TASK_DECL_4(int, final_measuring, QDD, bool*, BDDVAR*, BDDVAR);