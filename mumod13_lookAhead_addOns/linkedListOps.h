#ifndef _linkedListOps_h
#define _linkedListOps_h

typedef struct linkedlist{
  int *modeMotifCount;
  struct linkedlist *next;
}motifStruct;

typedef struct motifcont{
  motifStruct *motif;
}motifContainer;

motifStruct* createNode();
void freeNode(motifStruct*);
void freeMotif(motifStruct*);
void freeMotifs(motifContainer*, int);
motifStruct* initializeMotif(int);
motifContainer* initializeMotifs(int*, int);
motifStruct* addNodeFront(motifStruct*, motifStruct*);
motifStruct* addNodeEnd(motifStruct*, motifStruct*);
motifStruct* delNodeFront(motifStruct*);
motifStruct* delNodeEnd(motifStruct*);
motifStruct* getLastNode(motifStruct*);

#endif
