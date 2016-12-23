#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<limits.h>
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "motifops.h"


double motifAddLeftScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, t, n, *start;
  double *values, s;
  /* for(i = 0; i < ds->n; i++){ */
  /*   if(labels[i] == mode && startPos[i] <= 0) return 0; */
  /*   if((ds->data)[i][startPos[i] - 1] > 3) return 0; */
  /* } */
  t = (m->t)[mode];
  n = m->n;
  values = (double*)malloc(sizeof(double)*m->featureValues);
  if(!values) printMessages(0, NULL);
  start = (int*)malloc(sizeof(int)*ds->n);
  if(!start) printMessages(0, NULL);
  copyLabels(startPos, start, ds->n);
  for(i = 0; i < m->featureValues; i++) values[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] == 0){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    if((ds->data)[i][startPos[i] - 1] > 3){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    values[(ds->data)[i][startPos[i] - 1]]++;
  }
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(start[i] < 0) continue;
    if(background[i][(ds->data)[i][start[i] - 1]] < 0.00001) continue;
    s = s + log(values[(ds->data)[i][start[i] - 1]] + 0.5) - log(t + m->featureValues*0.5);
    s = s - log(background[i][(ds->data)[i][start[i] - 1]]);
  }
  free(values);
  free(start);
  //  printf("Left add: %lf\n", s);
  s = s - m->lambda;

  //  if(s > 700) s = 700;
  /* printf("Left add 1: %lf\n", s); */
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  /* else return s; */
  
  
  //  printf("Left add 2: %lf\n", s);
  //  if(s > DBL_MAX) s = DBL_MAX;
  return s;
}

double motifRemoveLeftScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i;
  double s;
  if((m->mWidth)[mode] <= 1) return 0;
  /* for(i = 0; i < ds->n; i++) if(labels[i] == mode && startPos[i] < 0) return 0; */
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(background[i][startPos[i]] < 0.00001) continue;
    s = s + log(background[i][startPos[i]]);
    s = s - log((((m->motifs)[mode].motif)->modeMotifCount)[(ds->data)[i][startPos[i]]] + 0.5) + log((m->t)[mode] + m->featureValues*0.5);
  }
  s = s + m->lambda;
  //  printf("S value::: %lf\n", s);
  if(s > 700) s = 700;
  /* printf("Left remove 1: %lf\n", s); */
  s = pow(exp(1), s);

  if(isinf(s)) return (double)INT_MAX;
  //  else return s;

  //  if(s > DBL_MAX) s = DBL_MAX;
  return s;
}

double motifAddRightScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, t, n, *start;
  double *values, s;

  /* for(i = 0; i < ds->n; i++){ */
  /*   if(labels[i] == mode && (startPos[i] == -1 || (startPos[i] + (m->mWidth)[mode] >= (ds->features)[i]))) return 0; */
  /*   if((ds->data)[i][startPos[i] + (m->mWidth)[mode]] > 3) return 0; */
  /* } */
  t = (m->t)[mode];
  n = m->n;
  values = (double*)malloc(sizeof(double)*m->featureValues);
  if(!values) printMessages(0, NULL);
  start = (int*)malloc(sizeof(int)*ds->n);
  if(!start) printMessages(0, NULL);
  copyLabels(startPos, start, ds->n);
  for(i = 0; i < m->featureValues; i++) values[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] + (m->mWidth)[mode] >= (ds->features)[i]){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    if((ds->data)[i][startPos[i] + (m->mWidth)[mode]] > 3){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    values[(ds->data)[i][startPos[i] + (m->mWidth)[mode]]]++;
  }
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(start[i] < 0) continue;
    if(background[i][start[i] + (m->mWidth)[mode]] < 0.00001) continue;
    s = s + log(values[(ds->data)[i][start[i] + (m->mWidth)[mode]]] + 0.5) - log(t + m->featureValues*0.5);
    s = s - log(background[i][start[i] + (m->mWidth)[mode]]);
  }
  free(values);
  free(start);
  s = s - m->lambda;
  //  if(s > 700) s = 700;
  /* printf("Right add 1: %lf\n", s); */
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  //  else return s;

  //  if(s > DBL_MAX) s = DBL_MAX;
  return s;
}

void checkMotifCount(dataSet *ds, model *m, int *labels, int *startPos, int mode, int width){
  int **motif, i, j, flag;
  motifStruct *mf;
  motif = (int**)malloc(sizeof(int*)*width);
  if(!motif) printMessages(0, NULL);
  for(i = 0; i < width; i++){
    motif[i] = (int*)malloc(sizeof(int)*m->featureValues);
    if(!motif[i]) printMessages(0, NULL);
    for(j = 0; j < m->featureValues; j++) motif[i][j] = 0;
  }
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    for(j = 0; j < width; j++)
      motif[j][(ds->data)[i][j + startPos[i]]]++;
  }
  mf = (m->motifs)[mode].motif;
  i = 0;
  while(mf != NULL){
    i++;
    mf = mf->next;
  }
  /* printf("Widths: %d, %d\n", i, width); */
  flag = 0;
  int i1, j1;
  mf = (m->motifs)[mode].motif;
  for(i = 0; i < width; i++){
    for(j = 0; j < m->featureValues; j++){
      /* printf("%d\t%d\n", (mf->modeMotifCount)[j], motif[i][j]); */
      if((mf->modeMotifCount)[j] != motif[i][j]){
	/* flag = 1; */
	printf("Not matching for i: %d, j: %d, mode: %d, mWidth: %d, width: %d, total: %d!!\n%d\t%d\n\n", i, j, mode, (m->mWidth)[mode], width, (m->t)[mode], (mf->modeMotifCount)[j], motif[i][j]);
	for(i1 = 0; i1 < ds->n; i1++){
	  if(labels[i1] != mode) continue;
	  if(startPos[i1] < 0) continue;
	  printf("Index: %d\tStart: %d\t", i1, startPos[i1]);
	  for(j1 = startPos[i1]; j1 < startPos[i1] + width; j1++)
	    printf("%d\t", (ds->data)[i1][j1]);
	  printf("\n");
	}
	exit(0);
      }
    }
    mf = mf->next;
  }
  /* if(flag == 1){ printf("Not matching!\n"); exit(0);} */
  for(i = 0; i < width; i++) free(motif[i]);
  free(motif);
}

double motifRemoveRightScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i;
  double s;
  //  int *values;
  motifStruct *mf;
  if((m->mWidth)[mode] <= 1) return 0;
  /* for(i = 0; i < ds->n; i++) if(labels[i] == mode && startPos[i] < 0) return 0; */
  s = 0;
  //  printf("HERE 0\n");
  mf = getLastNode((m->motifs)[mode].motif);
  /* mf = (m->motifs)[mode].motif; */
  /* while(mf->next != NULL){ */
  /*   printf("Mode count 0: %d\n", (mf->modeMotifCount)[0]); */
  /*   mf = mf->next; */
  /* } */
  /* printf("Mode count 1: %d\n", (mf->modeMotifCount)[1]); */
  //  printf("%d\n", (mf->modeMotifCount)[0]);
  /* if((m->motifs)[mode].motif == NULL) printf("HERE NULL 1!!\n"); */
  /* if(getLastNode((m->motifs)[mode].motif) == NULL) printf("HERE NULL!!\n"); */
  //  printf("HERE 0_2\n");
  //  values = mf->modeMotifCount;
  //  printf("HERE 0_1\n");
  /* printf("Values:\n"); */
  /* for(i = 0; i < m->featureValues; i++) printf("%d\t", values[i]); */
  /* printf("\n"); */
  /* if(mode == 3) printf("Inside count: %d, %d\n", (mf->modeMotifCount)[2], (m->t)[mode]); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    //    printf("Value: %d\n", startPos[i] + (m->mWidth)[mode] - 1);
    if(startPos[i] < 0) continue;
    if(background[i][startPos[i] + (m->mWidth)[mode] - 1] < 0.00001) continue;
    /* printf("C: %d\n", (ds->data)[i][startPos[i] + (m->mWidth)[mode] - 1]); */
    s = s + log(background[i][startPos[i] + (m->mWidth)[mode] - 1]);
    s = s - log((mf->modeMotifCount)[(ds->data)[i][startPos[i] + (m->mWidth)[mode] - 1]] + 0.5) + log((m->t)[mode] + m->featureValues*0.5);
  }
  //  printf("Values 1\n");
  s = s + m->lambda;
  if(s > 700) s = 700;
  /* printf("Right remove 1: %lf\n", s); */
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  //  else return s;

  //  if(s > DBL_MAX) s = DBL_MAX;
  return s;
}

double motifShiftLeftScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, t, n, *start;
  motifStruct *mf;
  double *values, s;
  mf = getLastNode((m->motifs)[mode].motif);
  t = (m->t)[mode];
  n = m->n;
  values = (double*)malloc(sizeof(double)*m->featureValues);
  if(!values) printMessages(0, NULL);
  start = (int*)malloc(sizeof(int)*ds->n);
  copyLabels(startPos, start, ds->n);
  for(i = 0; i < m->featureValues; i++) values[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] == 0){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    if((ds->data)[i][startPos[i] - 1] > 3){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    values[(ds->data)[i][startPos[i] - 1]]++;
  }
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(start[i] < 0) continue;
    if(background[i][(ds->data)[i][start[i] - 1]] < 0.00001) continue;
    s = s + log(values[(ds->data)[i][start[i] - 1]] + 0.5) - log(t + m->featureValues*0.5);
    s = s - log(background[i][(ds->data)[i][start[i] - 1]]);
    s = s + log(background[i][(ds->data)[i][start[i] + (m->mWidth)[mode] - 1]]);
    s = s - log((mf->modeMotifCount)[(ds->data)[i][start[i] + (m->mWidth)[mode] - 1]] + 0.5) + log((m->t)[mode] + m->featureValues*0.5);
  }
  free(values);
  free(start);
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  //  else return s;

  return s;
}

double motifShiftRightScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  motifStruct *mf;
  int i, t, n, *start;
  double *values, s;
  t = (m->t)[mode];
  n = m->n;
  mf = (m->motifs)[mode].motif;
  values = (double*)malloc(sizeof(double)*m->featureValues);
  if(!values) printMessages(0, NULL);
  start = (int*)malloc(sizeof(int)*ds->n);
  if(!start) printMessages(0, NULL);
  copyLabels(startPos, start, ds->n);
  for(i = 0; i < m->featureValues; i++) values[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] + (m->mWidth)[mode] >= (ds->features)[i]){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    if((ds->data)[i][startPos[i] + (m->mWidth)[mode]] > 3){
      start[i] = -1;
      t--;
      n--;
      continue;
    }
    values[(ds->data)[i][startPos[i] + (m->mWidth)[mode]]]++;
  }
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(start[i] < 0) continue;
    if(background[i][start[i] + (m->mWidth)[mode]] < 0.00001) continue;
    s = s + log(values[(ds->data)[i][start[i] + (m->mWidth)[mode]]] + 0.5) - log(t + m->featureValues*0.5);
    s = s - log(background[i][start[i] + (m->mWidth)[mode]]);
    s = s + log(background[i][start[i] + 1]);
    s = s - log((mf->modeMotifCount)[(ds->data)[i][start[i]]] + 0.5) + log((m->t)[mode] + 0.5*m->featureValues);
  }
  free(values);
  free(start);
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  //  else return s;

  return s;
}

void motifShiftLeft(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, j;
  motifStruct *mf, *mff;
  mf = createNode();
  for(i = 0; i < m->featureValues; i++) (mf->modeMotifCount)[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] == 0){
      mff = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mff->modeMotifCount)[(ds->data)[i][j]]--;
	mff = mff->next;
      }
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    if((ds->data)[i][startPos[i] - 1] > 3){
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    (mf->modeMotifCount)[(ds->data)[i][startPos[i] - 1]]++;
    startPos[i]--;
  }
  (m->motifs)[mode].motif = delNodeEnd((m->motifs)[mode].motif);
  (m->motifs)[mode].motif = addNodeFront(mf, (m->motifs)[mode].motif);
}

void motifShiftRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, j;
  motifStruct *me, *mf;
  me = createNode();
  for(i = 0; i < m->featureValues; i++) (me->modeMotifCount)[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] + (m->mWidth)[mode] >= (ds->features)[i]){
      mf = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mf->modeMotifCount)[(ds->data)[i][j]]--;
	mf = mf->next;
      }
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    if((ds->data)[i][startPos[i] + (m->mWidth)[mode]] > 3){
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    (me->modeMotifCount)[(ds->data)[i][startPos[i] + (m->mWidth)[mode]]]++;
    startPos[i]++;
  }
  (m->motifs)[mode].motif = delNodeFront((m->motifs)[mode].motif);
  (m->motifs)[mode].motif = addNodeEnd(me, (m->motifs)[mode].motif);
}

void motifLeftIncrease(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, j;
  motifStruct *mf, *mff;
  mf = createNode();
  for(i = 0; i < m->featureValues; i++) (mf->modeMotifCount)[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] == 0){
      
      mff = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mff->modeMotifCount)[(ds->data)[i][j]]--;
	mff = mff->next;
      }

      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    if((ds->data)[i][startPos[i] - 1] > 3){
      mff = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mff->modeMotifCount)[(ds->data)[i][j]]--;
	mff = mff->next;
      }
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    (mf->modeMotifCount)[(ds->data)[i][startPos[i] - 1]]++;
    startPos[i]--;
  }
  (m->motifs)[mode].motif = addNodeFront(mf, (m->motifs)[mode].motif);
  (m->mWidth)[mode]++;
}

void motifLeftDecrease(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    startPos[i]++;
  }
  (m->motifs)[mode].motif = delNodeFront((m->motifs)[mode].motif);
  (m->mWidth)[mode]--;
}

void motifRightIncrease(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, j;
  motifStruct *me, *mf;
  me = createNode();
  for(i = 0; i < m->featureValues; i++) (me->modeMotifCount)[i] = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(startPos[i] + (m->mWidth)[mode] >= (ds->features)[i]){

      mf = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mf->modeMotifCount)[(ds->data)[i][j]]--;
	mf = mf->next;
      }

      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    if((ds->data)[i][startPos[i] + (m->mWidth)[mode]] > 3){
      mf = (m->motifs)[mode].motif;
      for(j = startPos[i]; j < startPos[i] + (m->mWidth)[mode]; j++){
	(mf->modeMotifCount)[(ds->data)[i][j]]--;
	mf = mf->next;
      }
      startPos[i] = -1;
      (m->t)[mode]--;
      m->n--;
      continue;
    }
    (me->modeMotifCount)[(ds->data)[i][startPos[i] + (m->mWidth)[mode]]]++;
  }
  (m->motifs)[mode].motif = addNodeEnd(me, (m->motifs)[mode].motif);
  (m->mWidth)[mode]++;
}

void motifRightDecrease(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  (m->motifs)[mode].motif = delNodeEnd((m->motifs)[mode].motif);
  (m->mWidth)[mode]--;
}
