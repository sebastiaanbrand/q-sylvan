#include "QASM_to_Sylvan.h"

#include "sylvan.h"
#include "sylvan_qdd_complex.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// void remove_multiline_comments(char *str)
// {
//     char *c = str;
//     c=strchr(c,'/');
//     if (c != NULL) {
//         if (str[c-str+1] == '*') {

//         }
//     }
// }

QDD handle_tokens(QDD qdd, char **tokens)
{
    char *temp;
    if(strcmp(tokens[0], "qreg") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        qdd = qdd_create_all_zero_state(atoi(temp));
    }
    else if(strcmp(tokens[0], "h") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        qdd = qdd_gate(qdd, GATEID_H, atoi(temp));
    }
    else if(strcmp(tokens[0], "x") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        qdd = qdd_gate(qdd, GATEID_X, atoi(temp));
    }
    else if(strcmp(tokens[0], "y") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        qdd = qdd_gate(qdd, GATEID_Y, atoi(temp));
    }
    else if(strcmp(tokens[0], "z") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        qdd = qdd_gate(qdd, GATEID_Z, atoi(temp));
    }
    else if(strcmp(tokens[0], "cx") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        char *temp2 = strtok(NULL, "[");
        temp2 = strtok(NULL, "]");
        printf("arg2: %d\n",atoi(temp2));
        qdd = qdd_cgate(qdd, GATEID_X, atoi(temp), atoi(temp2));
    }
    else if(strcmp(tokens[0], "measure") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    printf("arg1: %d\n",atoi(temp));
    return qdd;

}

void read_QASM(char *filename)
{
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL)
    {
        perror("Error while opening the file.\n");
        exit(-1);
    }
    
    char *line = NULL, *c;
    size_t len = 0;
    ssize_t read;
    char *tokens[2];
    QDD qdd = qdd_create_all_zero_state(0);
    while ((read = getline(&line, &len, f)) != -1) {
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t')) line++;
        // remove empty lines, trailing information after ';'
        c=strchr(line,';');
        printf("line: %s",line);
        if (c != NULL) {
            line[c-line] = '\0';
            // tokenize string
            tokens[0] = strtok(line, " ");
            tokens[1] = strtok(NULL, "");
            qdd = handle_tokens(qdd, tokens);
        }
    }
    fclose(f);
}

int main(int argc, char *argv[])
{
    // check for file
    char *filename = "";
    if(argc == 2)
    {
        filename = argv[1];
    }
    else
    {
        printf("%s: Exactly two arguments is needed!\n", argv[0]);
        return -1;
    }

    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s)\n", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<19);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    read_QASM(filename);

    sylvan_quit();
    lace_exit();

    return 0;
}