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

void handle_tokens(char **tokens, int size)
{
    char *temp;
    if(strcmp(tokens[0], "qreg") == 0)
    {
        temp = strtok(tokens[1], "[]");
        temp = strtok(NULL, "[]");
        // QDD qdd = qdd_create_all_zero_state(atoi(temp));
    }
    // elif(strcmp(tokens[0], "h") == 0)
    // {
    //     printf("%s",tokens[0]);
    // }
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
    char *tokens[16];
    int i = 0, j = 0, k = 0;

    while ((read = getline(&line, &len, f)) != -1) {
        // remove leading spaces
        while ((*line == ' ') || (*line == '\t')) line++;
        // remove empty lines, trailing information after ';'
        c=strchr(line,';');
        if (c != NULL) {
            line[c-line] = '\0';
            // tokenize string
            i = 0;
            tokens[i] = strtok(line, " ");
            while(tokens[i]!=NULL)
            {
                tokens[++i] = strtok(NULL, " ");
            }
            handle_tokens(tokens, i);
            // printf("i: %d\n",i);
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