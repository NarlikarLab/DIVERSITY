#include<stdio.h>
#include<stdlib.h>

/* Print error messages and exit program */
void printMessages(int i, char *s){
  switch(i){
  case 0: printf("ERROR: Invalid malloc operation\n"); exit(1); break;
  case 1: printf("ERROR: Invalid fasta file %s\n", s); exit(1); break;
  case 2: printf("ERROR: 0 number of sequences in line %s\n", s); exit(1); break;
  case 3: printf("ERROR: Invalid fasta file %s. Cannot have both T's and U's\n", s); exit(1); break;
  case 4: printf("ERROR: Could not open file %s\n", s); exit(1); break;
  case 5: printf("ERROR: Not enough information present in FASTA file. Please check if any sequence contains too many N's\n"); exit(1); break;
  }
}
