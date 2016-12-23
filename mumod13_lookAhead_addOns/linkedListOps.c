#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "linkedListOps.h"
#include "messages.h"

motifStruct* createNode(){
  motifStruct *m = (motifStruct*)malloc(sizeof(motifStruct));
  if(!m) printMessages(0, NULL);
  m->modeMotifCount = (int*)malloc(sizeof(int)*4);
  if(!(m->modeMotifCount)) printMessages(0, NULL);
  m->next = NULL;
  return m;
}

void freeNode(motifStruct *m){
  free(m->modeMotifCount);
  m->next = NULL;
  m->modeMotifCount = NULL;
  free(m);
}

void freeMotif(motifStruct *m){
  motifStruct *mc;
  while(m != NULL){
    mc = m;
    m = m->next;
    freeNode(mc);
  }
}

void freeMotifs(motifContainer *m, int n){
  int i;
  for(i = 0; i < n; i++){
    freeMotif(m[i].motif);
    m[i].motif = NULL;
  }
  free(m);
}

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

motifContainer* initializeMotifs(int *widths, int n){
  motifContainer *mc;
  int i;
  mc = (motifContainer*)malloc(sizeof(motifContainer)*n);
  if(!mc) printMessages(0, NULL);
  for(i = 0; i < n; i++)
    mc[i].motif = initializeMotif(widths[i]);
  return mc;
}

motifStruct* addNodeFront(motifStruct *f, motifStruct *m){
  f->next = m;
  return f;
}

motifStruct* addNodeEnd(motifStruct *e, motifStruct *m){
  motifStruct *mc;
  mc = m;
  while(mc->next != NULL) mc = mc->next;
  mc->next = e;
  return m;
}

motifStruct* delNodeFront(motifStruct *m){
  motifStruct *mc;
  mc = m->next;
  freeNode(m);
  return mc;
}

motifStruct* delNodeEnd(motifStruct *m){
  motifStruct *mc, *mc1;
  mc = m;
  while((mc->next)->next != NULL) mc = mc->next;
  mc1 = mc->next;
  freeNode(mc1);
  mc->next = NULL;
  return(m);
}

motifStruct* getLastNode(motifStruct *m){
  motifStruct *mc;
  mc = m;
  //  printf("here!!\n");
  while(mc->next != NULL){
    //    printf("HERE000\n");
    mc = mc->next;
  }
  //  printf("here1!!\n");
  //  printf("Last node 0 count: %d\n", (mc->modeMotifCount)[0]);
  return mc;
}
