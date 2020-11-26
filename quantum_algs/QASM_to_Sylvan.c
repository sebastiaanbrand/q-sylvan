#include "sylvan_qdd_complex.h"

void read_QASM(string filename) {

   FILE *fp;
   fp = fopen(filename, "r");
   if (fp == NULL)
   {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
   }


}