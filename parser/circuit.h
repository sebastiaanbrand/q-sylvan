#include <stdbool.h>

#include "QASM_to_circuit.h"
#include "sylvan.h"

typedef struct Circuit {
    Gate** circuit;
    BDDVAR nvars;
    BDDVAR depth;
    BDDVAR* progress;
    QDD qdd;
} Circuit;

Circuit* create_circuit(char* filename);
void print_circuit(Circuit* c_s, bool vertical, bool show_rotation);
void delete_circuit(Circuit* c_s);
void skip_gate(Circuit* circuit_s, BDDVAR i);


#define advance(circuit_s,q) (CALL(advance,circuit_s,q));
TASK_DECL_2(bool, advance, Circuit*, BDDVAR);
#define get_gateid(gate) (CALL(get_gateid,gate));
TASK_DECL_1(BDDVAR, get_gateid, Gate);
