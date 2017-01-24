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
#include<float.h>
#include<limits.h>
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "motifops.h"

/* Compute score of increasing the motif width by one on the left */
double motifAddLeftScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, t, n, *start;
  double *values, s;
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
  s = s - m->lambda;
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  return s;
}

/* Compute score of decreasing the motif width by one on the left */
double motifRemoveLeftScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i;
  double s;
  if((m->mWidth)[mode] <= 1) return 0;
  s = 0;
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(background[i][startPos[i]] < 0.00001) continue;
    s = s + log(background[i][startPos[i]]);
    s = s - log((((m->motifs)[mode].motif)->modeMotifCount)[(ds->data)[i][startPos[i]]] + 0.5) + log((m->t)[mode] + m->featureValues*0.5);
  }
  s = s + m->lambda;
  if(s > 700) s = 700;
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  return s;
}

/* Compute score of increasing the motif width by one on the right */
double motifAddRightScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i, t, n, *start;
  double *values, s;
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
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  return s;
}


/* Compute score of decreasing the motif width by one on the right */
double motifRemoveRightScore(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  int i;
  double s;
  motifStruct *mf;
  if((m->mWidth)[mode] <= 1) return 0;
  s = 0;
  mf = getLastNode((m->motifs)[mode].motif);
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    if(startPos[i] < 0) continue;
    if(background[i][startPos[i] + (m->mWidth)[mode] - 1] < 0.00001) continue;
    s = s + log(background[i][startPos[i] + (m->mWidth)[mode] - 1]);
    s = s - log((mf->modeMotifCount)[(ds->data)[i][startPos[i] + (m->mWidth)[mode] - 1]] + 0.5) + log((m->t)[mode] + m->featureValues*0.5);
  }
  s = s + m->lambda;
  if(s > 700) s = 700;
  s = pow(exp(1), s);
  if(isinf(s)) return (double)INT_MAX;
  return s;
}

/* Increase motif size by one from left */
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

/* Decrease motif size by one from left */
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

/* Increase motif size by one from right */
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

/* Decrease motif size by one from right */
void motifRightDecrease(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode){
  (m->motifs)[mode].motif = delNodeEnd((m->motifs)[mode].motif);
  (m->mWidth)[mode]--;
}
