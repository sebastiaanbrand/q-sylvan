#include "QASM_to_Sylvan.h"
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

void handle_tokens(char **tokens)
{
    char *temp;
    if(strcmp(tokens[0], "qreg") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        // QDD qdd = qdd_create_all_zero_state(atoi(temp));
    }
    else if(strcmp(tokens[0], "h") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    else if(strcmp(tokens[0], "x") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    else if(strcmp(tokens[0], "y") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    else if(strcmp(tokens[0], "z") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    else if(strcmp(tokens[0], "cx") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
        char *temp2 = strtok(NULL, "[");
        temp2 = strtok(NULL, "]");
        printf("temp2: %d\n",atoi(temp2));
    }
    else if(strcmp(tokens[0], "measure") == 0)
    {
        printf("command: %s\n",tokens[0]);
        temp = strtok(tokens[1], "[");
        temp = strtok(NULL, "]");
    }
    else
    {
        return;
    }
    printf("temp1: %d\n",atoi(temp));

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
            printf("0: %s\n",tokens[0]);
            printf("1: %s\n",tokens[1]);
            handle_tokens(tokens);
            // for (j = 0; j < i; j++) {
            //     c = strchr(tokens[j], ';');
            //     printf("%d, %s\n",j, tokens[j]);
            // }
            // if( k == 4)
            //     break;
            // k++;
        }
    }
    fclose(f);
}