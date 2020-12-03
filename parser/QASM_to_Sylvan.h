#include <stdbool.h>

#include "sylvan.h"

void read_QASM(char *filename);

#define handle_intermediate_measure(qdd,measurements,nvars) (CALL(handle_intermediate_measure,qdd,measurements,nvars));
TASK_DECL_3(QDD, handle_intermediate_measure, QDD, bool*, BDDVAR);

#define handle_single_qubit_gate(qdd,target,gate_id) (CALL(handle_single_qubit_gate,qdd,target,gate_id));
TASK_DECL_3(QDD, handle_single_qubit_gate, QDD, char*, uint32_t);

#define handle_tokens(qdd,tokens,measurements,nvars) (CALL(handle_tokens,qdd,tokens,measurements,nvars));
TASK_DECL_4(QDD, handle_tokens, QDD, char**,bool*, BDDVAR*);

#define final_measuring(qdd,measurements,nvars) (CALL(final_measuring,qdd,measurements,nvars));
TASK_DECL_3(double*, final_measuring, QDD, bool*, BDDVAR*);