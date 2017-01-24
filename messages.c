///////////////////// DIVERSITY /////////////////////

//    DIVERSITY is a tool to explore multiple ways of protein-DNA
//    binding in the genome. More information can be found in the README file.
//    Copyright (C) 2015  Sneha Mitra, Anushua Biswas and Leelavati Narlikar

//    DIVERSITY is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    DIVERSITY is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////

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
