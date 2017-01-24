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
#include<math.h>
#include "linkedListOps.h"
#include "messages.h"

/* Creates a node to store counts at a position in the given motif. Counts stored in the form of a linked list */
motifStruct* createNode(){
  motifStruct *m = (motifStruct*)malloc(sizeof(motifStruct));
  if(!m) printMessages(0, NULL);
  m->modeMotifCount = (int*)malloc(sizeof(int)*4);
  if(!(m->modeMotifCount)) printMessages(0, NULL);
  m->next = NULL;
  return m;
}

/* Free a node of the linked list */
void freeNode(motifStruct *m){
  free(m->modeMotifCount);
  m->next = NULL;
  m->modeMotifCount = NULL;
  free(m);
}

/* Free the entire linked list storing count information of a given motif */
void freeMotif(motifStruct *m){
  motifStruct *mc;
  while(m != NULL){
    mc = m;
    m = m->next;
    freeNode(mc);
  }
}

/* Free all motifs that are stored in the form of a linked list */
void freeMotifs(motifContainer *m, int n){
  int i;
  for(i = 0; i < n; i++){
    freeMotif(m[i].motif);
    m[i].motif = NULL;
  }
  free(m);
}

/* Initialize a motif structure of given size */
motifStruct* initializeMotif(int size){
  int i;
  motifStruct *m, *mc;
  m = createNode();
  mc = m;
  for(i = 1; i < size; i++){
    mc->next = createNode();
    mc = mc->next;
  }
  return m;
}

/* Initialize al motifs for a given model */
motifContainer* initializeMotifs(int *widths, int n){
  motifContainer *mc;
  int i;
  mc = (motifContainer*)malloc(sizeof(motifContainer)*n);
  if(!mc) printMessages(0, NULL);
  for(i = 0; i < n; i++)
    mc[i].motif = initializeMotif(widths[i]);
  return mc;
}

/* Add a node to the front. Meaning motif width increased by one on the left */
motifStruct* addNodeFront(motifStruct *f, motifStruct *m){
  f->next = m;
  return f;
}

/* Add a node to the end. Meaning motif width increased by one on the right */
motifStruct* addNodeEnd(motifStruct *e, motifStruct *m){
  motifStruct *mc;
  mc = m;
  while(mc->next != NULL) mc = mc->next;
  mc->next = e;
  return m;
}

/* Delete a node from the front. Meaning motif width decreased by one on the left */
motifStruct* delNodeFront(motifStruct *m){
  motifStruct *mc;
  mc = m->next;
  freeNode(m);
  return mc;
}

/* Delete a node from the end. Meaning motif width decreased by one on the right */
motifStruct* delNodeEnd(motifStruct *m){
  motifStruct *mc, *mc1;
  mc = m;
  while((mc->next)->next != NULL) mc = mc->next;
  mc1 = mc->next;
  freeNode(mc1);
  mc->next = NULL;
  return(m);
}

/* Get pointer to the last node of the motif */
motifStruct* getLastNode(motifStruct *m){
  motifStruct *mc;
  mc = m;
  while(mc->next != NULL){
    mc = mc->next;
  }
  return mc;
}
