#include "circuit.h"

/**
 * Prints the QASM circuit given by <filename>.
 * 
 * PARAMETERS:
 * - filename: the path to the file containing the QASM circuit code
 * 
 * FLAGS:
 * [-o optimize] (optional) optimize the circuit if true. This option will remove negating gates
 */
int main(int argc, char *argv[])
{
    // Initialise flag parameters
    char *filename = argv[1];
    bool optimize = false;
    bool rotation = false;
    int opt;

    // Read flags from cmd and set parameters
    while((opt = getopt(argc, argv, "or")) != -1) {
        switch(opt) {
            case 'o':
                optimize = true;
                break;
            case 'r':
                rotation = true;
                break;
            default:
                fprintf(stderr, "usage: %s file [-o optimize][-r show_rotations]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    // Check if a file is given, if not, return an error
    if(access(filename, F_OK) != 0)
    {
        fprintf(stderr, "Invalid QASM file.\n");
        exit(EXIT_FAILURE);
    }

    // Create circuit struct
    C_struct c_s = make_c_struct(filename, optimize);
    // Print circuit
    print_c_struct(c_s, rotation);
    // Free circuit variables
    delete_c_struct(&c_s);
    return 0;
}