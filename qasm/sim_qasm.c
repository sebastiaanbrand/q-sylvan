#include "qsylvan.h"
#include "simple_parser.h"


int main(int argc, char *argv[])
{
    quantum_op_t* circuit = parse_qasm_file(argv[1]);
    
    // TODO: run circuit

    free_quantum_ops(circuit);
}
