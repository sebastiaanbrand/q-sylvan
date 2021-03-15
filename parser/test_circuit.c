#include "circuit.h"

int main(int argc, char *argv[])
{
    // check for file
    char *filename = "";
    bool optimize = false;
    int opt;
    while((opt = getopt(argc, argv, "f:o")) != -1) {
        switch(opt) {
            case 'f':
                filename = optarg;
                break;
            case 'o':
                optimize = true;
                break;
        }
    }
    if(strcmp(filename, "") == 0)
    {
        printf("Give filename of qasm file.\n");
        return 1;
    }

    C_struct c_s1 = make_c_struct(filename, optimize);
    // Circuit *c_s2 = create_circuit(filename);
    print_c_struct(c_s1, false, true);
    printf("\n\n");
    // print_circuit(c_s2, false, true);
    delete_c_struct(&c_s1);
    // delete_circuit(c_s2);
    return 0;
}