#include "qsylvan.h"
#include "simple_parser.h"


int main(int argc, char *argv[])
{
    quantum_circuit_t* circuit = parse_qasm_file(argv[1]);
    
    // TODO: run circuit
    print_quantum_circuit(circuit);
    free_quantum_circuit(circuit);
}
