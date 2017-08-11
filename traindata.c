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
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "motifops.h"
#include "traindata.h"

FILE *fp;
int citer;

/* Compute slope for given set likelihood values */
double linearSlope(queueStruct *q, double xsum, double ysum, double xxsum, double xysum, int n){
  double a1;
  a1 = (n * xysum - xsum * ysum) / (n * xxsum - xsum * xsum);
  return a1;
}

/* Enqueue a likelihood value */
void enqueue(queueStruct *queue, int *front, int *rear, int *isFull, int size, double value, double *sum, double *sumxy){
  if(*isFull == size) dequeue(queue, front, rear, size, sum, sumxy);
  if(*front == -1) *front = *rear = 0;
  else *rear = (*rear + 1)%size;
  queue[*rear].y = value;
  *sum = *sum + value;
  if(sumxy != NULL) *sumxy = *sumxy + value * (*rear + 1);
  if(*isFull < size) *isFull = *isFull + 1;
}

/* Dequeue a likelihood value */
void dequeue(queueStruct *queue, int *front, int *rear, int size, double *sum, double *sumxy){
  *sum = *sum - queue[*front].y;
  if(sumxy != NULL) *sumxy = *sumxy - queue[*front].y * (*front + 1);
  if(*front == *rear) *front = *rear = -1;
  else *front = (*front + 1)%size;
}

/* Remove sequence from data set */
void addRemoveDataPoint(model *m, dataSet *ds, int *labels, int *startPos, int index, int ar){
  int i, mode;
  motifStruct *m1;
  mode = labels[index];
  if(startPos[index] == -1) return;
  (m->t)[mode] = (m->t)[mode] + ar;
  m->n = m->n + ar;
  m1 = (m->motifs)[mode].motif;
  for(i = 0; i < (m->mWidth)[mode]; i++){
    if((i + startPos[index]) >= (ds->features)[index]) break;
    (m1->modeMotifCount)[(ds->data)[index][i + startPos[index]]] = (m1->modeMotifCount)[(ds->data)[index][i + startPos[index]]] + ar;
    m1 = m1->next;
  }
}

/* Compute parameters that maximize the posterior */
int bestPosteriorParameters(dataSet *ds, model *m, double *background, int *labels, int *startPos, int index, int *t, int n, int fv){
  int i, j, c, id;
  double s, s1, max;
  max = -1;
  id = 0;
  c = 0;
  for(i = 0; i < m->mode; i++){
    s1 = (double)(t[i] + m->alpha)/(n + m->mode*m->alpha);
    for(j = 0; j < (ds->features)[index]; j++){
      s = likelihoodXi(m, ds, i, j, index, background);
      if(s < 0.0000000000001) s = -1;
      else s = s * s1;
      if(s > max){
    	max = s;
    	id = c;
      }
      c++;
    }
  }
  if(m->zoops == 1){
    if(max < 1) id = -1;
  }
  else if(m->zoops > 0){
    if((max * (1 - m->zoops) / (ds->features)[index]) < m->zoops) id = -1;
  }
  if(max == -1) id = -1;
  return id;
}

/* Update the best model */
int updateBestModel(dataSet *ds, model *m, double **background, int *labels, int *startPos, int fv){
  int i, c, flag, start, label, n;
  int *tmpCount;
  tmpCount = (int*)malloc(sizeof(int)*m->mode);
  if(!tmpCount) printMessages(0, NULL);
  for(i = 0; i < m->mode; i++)
    tmpCount[i] = (m->t)[i];
  n = m->n;
  flag = 0;
  for(i = 0; i < ds->n; i++){
    label = labels[i];
    start = startPos[i];
    c = bestPosteriorParameters(ds, m, background[i], labels, startPos, i, tmpCount, n, fv);
    if(c == -1){
      start = c;
      label = c;
    }
    else{
      label = c / (ds->features)[i];
      start = c % (ds->features)[i];
    }
    if((label != labels[i] || start != startPos[i]) /* && tmpCount[labels[i]] > 20 */){
      if(label != labels[i] && label == -1) n--;
      else if(label != labels[i] && labels[i] == -1) n++;
      if(labels[i] >= 0 && labels[i] < m->mode) tmpCount[labels[i]]--;
      labels[i] = label;
      if(labels[i] >= 0 && labels[i] < m->mode) tmpCount[labels[i]]++;
      startPos[i] = start;
      flag = 1;
    }
  }
  free(tmpCount);
  return flag;
}

/* Compute posterior */
double calculateFullPosterior(dataSet *ds, model *m, double **background, int *labels, int *startPos, int **featureCounts){
  double s, s1;
  int i, j, k;
  motifStruct *m1;
  s = 0;
  s1 = 0;
  for(i = 0; i < m->mode; i++){
    s = s + ((m->t)[i])*(log((m->t)[i] + m->alpha) - log(m->n + m->mode*m->alpha));
    s = s - m->lambda*(m->mWidth)[i];
    m1 = (m->motifs)[i].motif;
    for(j = 0; j < (m->mWidth)[i]; j++){
      for(k = 0; k < m->featureValues; k++){
	s1 = s1 + log((m1->modeMotifCount)[k] + 0.5) - log((m->t)[i] + m->featureValues*0.5);
      }
      m1 = m1->next;
    }
  }
  for(i = 0; i < ds->n; i++){
    if(startPos[i] == -1){
      for(j = 0; j < (ds->features)[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]);
      if(m->zoops > 0 && m->zoops < 1) s = s + log(m->zoops);
      else s = s - log((ds->features)[i]);
      continue;
    }
    if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][labels[i]]);
    else if(m->zoops == 1) s = s - log(featureCounts[i][labels[i]] + 1);
    else s = s - log(featureCounts[i][labels[i]]);
    m1 = (m->motifs)[labels[i]].motif;
    for(j = 0; j < (m->mWidth)[labels[i]]; j++){
      s = s + log((m1->modeMotifCount)[(ds->data)[i][j + startPos[i]]] + 0.5) - log((m->t)[labels[i]] + m->featureValues*0.5);
      m1 = m1->next;
    }
    for(j = 0; j < startPos[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]);
    for(j = startPos[i] + (m->mWidth)[labels[i]]; j < (ds->features)[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]);
  }
  s1 = s1 * (0.5 - 1);
  s = s + s1;
  return s;
}

/* Calculate likelihood */
double calculateLikelihood(model *m, dataSet *ds, int *lp, int *startPos, double **background, int **featureCounts){
  double s, s1, lxi;
  int i, j, k;
  motifStruct *m1;
  s = 0;
  s1 = 0;

  for(i = 0; i < m->mode; i++){
    s = s + ((m->t)[i] + m->alpha)*(log((m->t)[i] + m->alpha) - log(m->n + m->mode*m->alpha));
    s = s - m->lambda*(m->mWidth)[i];
  }
  for(i = 0; i < ds->n; i++){
    lxi = 0;
    if(startPos[i] != -1) lxi = likelihoodXiFaster(m, ds, lp[i], startPos[i], i, background[i]);
    if(lxi != 0){
      lxi = log(lxi);
      if(m->zoops == 0) lxi = lxi - log(featureCounts[i][lp[i]]);
      else if(m->zoops == 1) lxi = lxi - log(featureCounts[i][lp[i]] + 1);
      else lxi = lxi + log(1 - m->zoops) - log(featureCounts[i][lp[i]]);
      s = s + lxi;
    }
    else if(startPos[i] == -1 && m->zoops != 0) s = s + log(m->zoops);
    else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][lp[i]]);
  }
  for(i = 0; i < m->mode; i++){
    m1 = (m->motifs)[i].motif;
    for(j = 0; j < (m->mWidth)[i]; j++){
      for(k = 0; k < m->featureValues; k++){
	s1 = s1 + log((m1->modeMotifCount)[k] + 0.5) - log((m->t)[i] + m->featureValues*0.5);
      }
      m1 = m1->next;
    }
  }
  s1 = s1 * (0.5 - 1);
  s = s + s1;
  return s;
}

/* Calculate likelihood of all sequences belonging to the given mode */
double calculateLikelihoodMode(model *m, dataSet *ds, int *labels, int *startPos, double **background, int **featureCounts, int mode){
  double s, s1;
  int i, j;
  motifStruct *mf;
  s = ((m->t)[mode] + m->alpha)*(log((m->t)[mode] + m->alpha) - log(m->n + m->mode*m->alpha));
  s = s - m->lambda*(m->mWidth)[mode];
  for(i = 0; i < ds->n; i++){
    if(labels[i] != mode) continue;
    s1 = 0;
    if(startPos[i] != -1) s1 = likelihoodXiFaster(m, ds, mode, startPos[i], i, background[i]);
    if(s1 != 0){
      s1 = log(s1);
      if(m->zoops == 0) s1 = s1 - log(featureCounts[i][mode]);
      else if(m->zoops == 1) s1 = s1 - log(featureCounts[i][mode] + 1);
      else s1 = s1 + log(1 - m->zoops) - log(featureCounts[i][mode]);
      s = s + s1;
    }
    else if(startPos[i] == -1 && m->zoops != 0) s = s + log(m->zoops);
    else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][mode]);
  }
  mf = (m->motifs)[mode].motif;
  s1 = 0;
  for(i = 0; i < (m->mWidth)[mode]; i++){
    for(j = 0; j < m->featureValues; j++){
      s1 = s1 + log((mf->modeMotifCount)[j] + 0.5) - log((m->t)[mode] + m->featureValues*0.5);
    }
    mf = mf->next;
  }

  s1 = s1 * (0.5 - 1);
  s = s + s1;
  return s;
}

/* Calculate likelihood of a single sequence */
double likelihoodXi(model *m, dataSet *ds, int mode, int start, int index, double *background){
  int i;
  double s;
  motifStruct *m1;
  s = 1;
  if(start == -1) return 0;
  if((start + (m->mWidth)[mode]) > (ds->features)[index]) return 0;
  m1 = (m->motifs)[mode].motif;
  for(i = 0; i < (m->mWidth)[mode]; i++){
    if((ds->data)[index][i + start] > 3){
      return 0;
    }
    s = s * ((m1->modeMotifCount)[(ds->data)[index][i + start]] + 0.5) / ((m->t)[mode] + m->featureValues*0.5);
    m1 = m1->next;
    if(background[i + start] <= 0.00001) continue;
    s = s / (background[i + start]);
  }
  return s;
}

/* Calculate likelihood of a single sequence faster */
double likelihoodXiFaster(model *m, dataSet *ds, int mode, int start, int index, double *background){
  int i;
  double s;
  motifStruct *m1;
  s = 1;
  m1 = (m->motifs)[mode].motif;
  for(i = 0; i < (m->mWidth)[mode]; i++){
    s = s * ((m1->modeMotifCount)[(ds->data)[index][i + start]] + 0.5) / ((m->t)[mode] + m->featureValues*0.5);
    m1 = m1->next;
    if(background[i + start] <= 0.00001) continue;
    s = s / (background[i + start]);
  }
  return s;
}

/* Sample label for given sequence */
int sampleNewLabel(model *m, dataSet *ds, int *startPos, int index, unsigned int *seed, int flag, double *background, int *featureCounts){
  double *values;
  int i, j;
  values = (double*)malloc(sizeof(double)*m->mode);
  if(!values) printMessages(0, NULL);
  for(i = 0; i < m->mode; i++){
    values[i] = ((m->t)[i] + m->alpha) / (m->n + m->mode*m->alpha);
    values[i] = values[i] * likelihoodXi(m, ds, i, startPos[index], index, background);
    if(values[i] != 0){
      if(m->zoops == 0) values[i] = values[i] / featureCounts[i];
      else if(m->zoops == 1) values[i] = values[i] / (featureCounts[i] + 1);
      else values[i] = values[i] * (1 - m->zoops) / featureCounts[i];
    }
  }
  if(flag == 1) j = arrayMax(values, m->mode);
  else{
    j = sample(values, m->mode, ((double)rand_r(seed))/(RAND_MAX));
  }
  free(values);
  return j;
}

/* Sample start position of motif for given sequence */
int sampleStartPos(model *m, dataSet *ds, int mode, int index, unsigned int *seed, int flag, double *background){
  double *values, sum;
  int i, j, count;
  sum = 0;
  count = 0;
  values = (double*)malloc(sizeof(double)*((ds->features)[index] + (m->zoops > 0)));
  if(!values) printMessages(0, NULL);
  for(i = 0; i < (ds->features)[index] + (m->zoops > 0); i++) values[i] = 0;
  i = 0;
  while(i < (ds->features)[index]){
    if((ds->lookAhead)[index][i] < (m->mWidth)[mode]){
      i = i + (ds->lookAhead)[index][i] + 1;
      continue;
    }
    values[i] = likelihoodXiFaster(m, ds, mode, i, index, background);
    sum = sum + values[i];
    i++;
  }
  if(sum < 0.00001){
    free(values);
    return -1;
  }

  if(m->zoops > 0){
    if(m->zoops == 1) values[(ds->features)[index]] = 1;
    else if (sum == 0) values[(ds->features)[index]] = 1;
    else values[(ds->features)[index]] = sum * m->zoops / (1 - m->zoops);
  }
  if(count == (ds->features)[index]) j = -1;
  else if(flag == 1) j = arrayMax(values, (ds->features)[index]);
  else{
    j = sample(values, (ds->features)[index] + (m->zoops > 0), ((double)rand_r(seed))/(RAND_MAX));
  }
  free(values);
  if(j == (ds->features)[index]) j = -1;
  return j;
}

/* Compute the width increase/decrease of a given motif from the left */
int sampleMotifWidthLeft(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int minWidth, int maxWidth, unsigned int *seed){
  int i;
  double *values;

  values = (double*)malloc(sizeof(double)*3);
  if(!values) printMessages(0, NULL);
  values[0] = 1;
  values[1] = motifRemoveLeftScore(ds, m, background, labels, startPos, mode);
  values[2] = motifAddLeftScore(ds, m, background, labels, startPos, mode);
  if((m->mWidth)[mode] == minWidth) values[1] = 0;
  if(maxWidth != 0 && (m->mWidth)[mode] >= maxWidth) values[2] = 0;
  i = sample(values, 3, ((double)rand_r(seed))/(RAND_MAX));
  if(i == 1) motifLeftDecrease(ds, m, background, labels, startPos, mode);
  else if(i == 2) motifLeftIncrease(ds, m, background, labels, startPos, mode);

  free(values);
  return i;
}

/* Compute the width increase/decrease of a given motif from the right */
int sampleMotifWidthRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int minWidth, int maxWidth, unsigned int *seed){
  int i;
  double *values;
  values = (double*)malloc(sizeof(double)*3);
  if(!values) printMessages(0, NULL);
  values[0] = 1;
  values[1] = motifRemoveRightScore(ds, m, background, labels, startPos, mode);
  values[2] = motifAddRightScore(ds, m, background, labels, startPos, mode);  
  if((m->mWidth)[mode] == minWidth) values[1] = 0;
  if(maxWidth != 0 && (m->mWidth)[mode] >= maxWidth) values[2] = 0;
  i = sample(values, 3, ((double)rand_r(seed))/(RAND_MAX));
  if(i == 2) motifRightIncrease(ds, m, background, labels, startPos, mode);
  else if(i == 1) motifRightDecrease(ds, m, background, labels, startPos, mode);
  free(values);
  return i;
}

double EMLike(model *m, dataSet *ds, int *labels, int *startPos, double **background){
  int j;
  double maxLikelihood;
  model *mc;
  int *lpc, *spc;
  j = 0;
  lpc = (int*)malloc(sizeof(int)*ds->n);
  spc = (int*)malloc(sizeof(int)*ds->n);
  copyLabels(labels, lpc, ds->n);
  copyLabels(startPos, spc, ds->n);
  mc = createModel(m->mode, ds, lpc, spc, m->alpha, m->lambda, m->zoops, m->mWidth);
  maxLikelihood = calculateFullPosterior(ds, mc, background, lpc, spc, NULL);
  while(updateBestModel(ds, mc, background, lpc, spc, j) && j < 25){
    getMotifCount(mc, ds, spc, lpc);
    maxLikelihood = calculateFullPosterior(ds, mc, background, lpc, spc, NULL);
    j++;
  }
  freeModel(mc);
  free(lpc);
  free(spc);
  return maxLikelihood;
}

/* model training */
trainOut* trainData(dataSet *ds, int mode, float fast, float alpha, float lambda, double zoops, unsigned int seed, double **background, int *mWidth, int minWidth, int maxWidth, char *filename, char *likelihoodInfoFile){
  double maxLikelihood, tmpLikelihood, *modeLikes;
  int i, j, j1, k, oldLabel, oldStart, flag, count, iterations;
  int *lpc, *spc;
  motifStruct *m1;
  int *widths;
  int **labelcounts;
  double c1 = 1000, *x;
  int oldChangeIndex, motifWidthChange;
  int front, rear, qSize, movingCounter, movingOffset, isFull;
  double qSum, xsum, xxsum, xysum;
  queueStruct *queue;
  trainOut *to;
  int *labels, *startPos;
  int **featureCounts, **featureCountsCopy;
  model *m;
  FILE *fp3, *fp2;

  int i1, j2, k1;
  
  /* Initializations */
  featureCounts = getFeatureCounts(ds, mode, mWidth[0]);
  featureCountsCopy = (int**)malloc(sizeof(int*)*ds->n);
  if(!featureCountsCopy) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++)
    {
      featureCountsCopy[i] = (int*)malloc(sizeof(int)*mode);
      if(!featureCountsCopy[i]) printMessages(0, NULL);
      for(j = 0; j < mode; j++) featureCountsCopy[i][j] = featureCounts[i][j];
      
    }
  
  to = (trainOut*)malloc(sizeof(trainOut));
  if(!to) printMessages(0, NULL);
  labels = (int*)malloc(sizeof(int)*ds->n);
  if(!labels) printMessages(0, NULL);
  startPos = (int*)malloc(sizeof(int)*ds->n);
  if(!startPos) printMessages(0, NULL);

  initializeLabelStartPos(ds, labels, startPos, mode, mWidth, &seed);
  m = createModel(mode, ds, labels, startPos, alpha, lambda, zoops, mWidth);

  oldChangeIndex = -1;
  motifWidthChange = 0;
  qSize = ds->n;
  movingCounter = 0;
  movingOffset = ds->n/4;
  front = -1;
  rear = -1;
  xsum = (double) qSize * (qSize + 1) / 2;
  xxsum = (double) qSize * (qSize + 1) * (2 * qSize + 1) / 6;
  qSum = 0;
  xysum = 0;
  isFull = 0;
  queue = (queueStruct*)malloc(sizeof(queueStruct)*qSize);
  if(!queue) printMessages(0, NULL);


  x = (double*)malloc(sizeof(double)*qSize);
  for(i = 0; i < qSize; i++) x[i] = i + 1;
  
  /* fp2 = fopen("likecalc.txt", "w"); */
  

  fp = NULL;
  if(filename[0] != '0') fp = fopen(filename, "w");
  lpc = (int*)malloc(sizeof(int)*ds->n);
  if(!lpc) printMessages(0, NULL);
  spc = (int*)malloc(sizeof(int)*ds->n);
  if(!spc) printMessages(0, NULL);
  widths = (int*)malloc(sizeof(int)*m->mode);
  if(!widths) printMessages(0, NULL);

  labelcounts = (int**)malloc(sizeof(int*)*ds->n);
  if(!labelcounts) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++){
    labelcounts[i] = (int*)malloc(sizeof(int)*m->mode);
    if(!labelcounts[i]) printMessages(0, NULL);
    for(j = 0; j < m->mode; j++) labelcounts[i][j] = 0;
    labelcounts[i][labels[i]]++;
  }

  copyLabels(labels, lpc, ds->n);
  copyLabels(startPos, spc, ds->n);
  copyLabels(m->mWidth, widths, m->mode);

  modeLikes = (double*)malloc(sizeof(double)*m->mode);
  if(!modeLikes) printMessages(0, NULL);

  for(i = 0; i < m->mode; i++) modeLikes[i] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, i);
  j = 0;
  i = 0;
  maxLikelihood = calculateLikelihood(m, ds, lpc, spc, background, featureCounts);

  tmpLikelihood = maxLikelihood;

  enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);

  count = 0;
  flag = 0;
  iterations = ds->n/fast;
  //  if((fast >= 0 && fast < 0.00001) || (fast < 0 && fast > -0.00001)) iterations = 0;
  if((fast >= 0 && fast < 0.00001)) iterations = 0;
  
  /* Iterate over the entire dataset n times where n is the size of the data set */
  while(1){
    j++;
    if(j == iterations) break;
    /* Iterate over all sequences of data set */
    for(i = 0; i < ds->n; i++){
      oldLabel = lpc[i];
      oldStart = spc[i];
      addRemoveDataPoint(m, ds, lpc, spc, i, -1);

      /* Sample label and motif start posiiton for given sequence */
      if(lpc[i] != -1) spc[i] = sampleStartPos(m, ds, lpc[i], i, &seed, 0, background[i]);
      if(spc[i] != -1){
      	lpc[i] = sampleNewLabel(m, ds, spc, i, &seed, 0, background[i], featureCounts[i]);
	labelcounts[i][lpc[i]]++;
      }

      addRemoveDataPoint(m, ds, lpc, spc, i, 1);

      if(oldLabel == lpc[i] && oldStart == spc[i]){

	/* Check of all  likelihood values lie in a straight line to see if it has converged */
	enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
	movingCounter = (movingCounter + 1)%(movingOffset + 1);
	if(movingCounter == movingOffset && isFull == qSize){
	  movingCounter = 0;
	  c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	  if((c1 <= 0 && c1 > -0.001) || (c1 >= 0 && c1 < 0.001)){
	    flag = 1;
	  }
	}
	continue;
      }

      modeLikes[oldLabel] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, oldLabel);
      if(oldLabel != lpc[i]) modeLikes[lpc[i]] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, lpc[i]);
      tmpLikelihood = 0;
      for(k = 0; k < m->mode; k++) tmpLikelihood = tmpLikelihood + modeLikes[k];
      if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
      /* fprintf(fp2, "%lf\n", calculateLikelihood(m, ds, lpc, spc, background, featureCounts)); */
      /* Check of all  likelihood values lie in a straight line to see if it has converged */
      enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
      movingCounter = (movingCounter + 1)%(movingOffset + 1);
      if(movingCounter == movingOffset && isFull == qSize){
	movingCounter = 0;
	c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	if((c1 <= 0 && c1 > -0.001) || (c1 >= 0 && c1 < 0.001)){
	  flag = 1;
	}
      }

      /* Update best model if the current likelihood is the maximum likelihood */
      if(tmpLikelihood > maxLikelihood){
	/* printf("max: %lf, calc: %lf\n", tmpLikelihood, calculateLikelihood(m, ds, lpc, spc, background, featureCounts)); */
	
	copyLabels(lpc, labels, ds->n);
	copyLabels(spc, startPos, ds->n);
	copyLabels(m->mWidth, widths, m->mode);

	for(i1 = 0; i1 < ds->n; i1++)
	  for(k1 = 0; k1 < mode; k1++)
	    featureCountsCopy[i1][k1] = featureCounts[i1][k1];
	
	maxLikelihood = tmpLikelihood;

	  /* fp3 = fopen("learnedlabels.txt", "w"); */
	  /* for(k = 0; k < ds->n; k++) fprintf(fp3, "%d\n", labels[k]); */
	  /* fclose(fp3); */

	  /* fp3 = fopen("learnedstart.txt", "w"); */
	  /* for(k = 0; k < ds->n; k++) fprintf(fp3, "%d\n", startPos[k]); */
	  /* fclose(fp3); */

	  /* fp3 = fopen("learnedwidth.txt", "w"); */
	  /* for(k = 0; k < m->mode; k++) fprintf(fp3, "%d\n", widths[k]); */
	  /* fclose(fp3); */

	  /* printf("feature count:\n"); */
	  /* for(i1 = 0; i1 < m->mode; i1++) */
	  /*   { */
	  /*     m1 = (m->motifs)[i1].motif; */
	  /*     for(j2 = 0; j2 < (m->mWidth)[i1]; j2++) */
	  /* 	{ */
	  /* 	  printf("%d\t", (m1->modeMotifCount)[1]); */
	  /* 	  m1 = m1->next; */
	  /* 	} */
	  /*     printf("\n"); */
	  /*   } */
	  
      }	
    }


    count++;
    if(j%(ds->n/10) != 0) continue;

    /* Update motif size after every n/10 iteration. So, the motif size would be changed at most 10 times while training */
    if(j < ds->n/10){
      oldChangeIndex = -1;
      motifWidthChange = 0;
    }
      /* Check of all  likelihood values lie in a straight line to see if it has converged */
    if(oldChangeIndex > 0 && motifWidthChange >= 2 && ((c1 <= 0 && c1 > -0.001) || (c1 >= 0 && c1 < 0.001))) break;

    oldChangeIndex = j;
    count = 0;

    /* Update motif size for every mode */
    for(i = 0; i < m->mode; i++){
      while(1){
	flag = 0;
	if((m->t)[i] < 10){ 
	  break;
	}
	else{
	  k = sampleMotifWidthLeft(ds, m, background, lpc, spc, i, minWidth, maxWidth, &seed);
	  if(k == 2) flag = 1;
	  j1 = sampleMotifWidthRight(ds, m, background, lpc, spc, i, minWidth, maxWidth, &seed);
	  if(j1 == 2) flag = 1;
	  k = k + j1;
	}
	if(k == 0){
	  motifWidthChange++;
	  break;
	}
	motifWidthChange = 0;
	modeFeatureCount(ds, featureCounts, (m->mWidth)[i], i);
	tmpLikelihood = calculateLikelihood(m, ds, lpc, spc, background, featureCounts);
	
	enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
	movingCounter = (movingCounter + 1)%(movingOffset + 1);
	if(movingCounter == movingOffset && isFull == qSize){
	  movingCounter = 0;
	  c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	}
	modeLikes[i] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, i);
	if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
	/* fprintf(fp2, "%lf\n", calculateLikelihood(m, ds, lpc, spc, background, featureCounts)); */

	/* Update best model if the current likelihood is the maximum likelihood */
	if(tmpLikelihood > maxLikelihood){
	  /* printf("max: %lf, calc: %lf\n", tmpLikelihood, calculateLikelihood(m, ds, lpc, spc, background, featureCounts)); */
	  maxLikelihood = tmpLikelihood;
	  copyLabels(lpc, labels, ds->n);
	  copyLabels(spc, startPos, ds->n);
	  copyLabels(m->mWidth, widths, m->mode);

	  for(i1 = 0; i1 < ds->n; i1++)
	    for(k1 = 0; k1 < mode; k1++)
	      featureCountsCopy[i1][k1] = featureCounts[i1][k1];
	  
	  /* fp3 = fopen("learnedlabels.txt", "w"); */
	  /* for(k = 0; k < ds->n; k++) fprintf(fp3, "%d\n", labels[k]); */
	  /* fclose(fp3); */

	  /* fp3 = fopen("learnedstart.txt", "w"); */
	  /* for(k = 0; k < ds->n; k++) fprintf(fp3, "%d\n", startPos[k]); */
	  /* fclose(fp3); */

	  /* fp3 = fopen("learnedwidth.txt", "w"); */
	  /* for(k = 0; k < m->mode; k++) fprintf(fp3, "%d\n", widths[k]); */
	  /* fclose(fp3); */

	  /* printf("motif count:\n"); */
	  /* for(i1 = 0; i1 < m->mode; i1++) */
	  /*   { */
	  /*     m1 = (m->motifs)[i1].motif; */
	  /*     for(j2 = 0; j2 < (m->mWidth)[i1]; j2++) */
	  /* 	{ */
	  /* 	  printf("%d\t", (m1->modeMotifCount)[1]); */
	  /* 	  m1 = m1->next; */
	  /* 	} */
	  /*     printf("\n"); */
	  /*   } */
	  
	}		
	if(flag == 0) break;
      }
    }
  }
  

  copyLabels(widths, m->mWidth, m->mode);
  freeMotifs(m->motifs, m->mode);
  m->motifs = initializeMotifs(m->mWidth, m->mode);
  getMotifCount(m, ds, startPos, labels);

  for(i1 = 0; i1 < ds->n; i1++)
    for(k1 = 0; k1 < mode; k1++)
      featureCounts[i1][k1] = featureCountsCopy[i1][k1];

  
  /* printf("Final motif count:\n"); */
  /* for(i1 = 0; i1 < m->mode; i1++) */
  /*   { */
  /*     m1 = (m->motifs)[i1].motif; */
  /*     for(j1 = 0; j1 < (m->mWidth)[i1]; j1++) */
  /* 	{ */
  /* 	  printf("%d\t", (m1->modeMotifCount)[1]); */
  /* 	  m1 = m1->next; */
  /* 	} */
  /*     printf("\n"); */
  /*   } */
  
  
  //////////////////////////////////////////////
  /* printf("OUTSIDE: %lf\n", maxLikelihood); */
  /* printf("Calculated: %lf\n", calculateLikelihood(m, ds, labels, startPos, background, featureCounts)); */
  
  /* fp3 = fopen("labels.txt", "w"); */
  /* for(i = 0; i < ds->n; i++) fprintf(fp3, "%d\n", labels[i]); */
  /* fclose(fp3); */

  /* fp3 = fopen("start.txt", "w"); */
  /* for(i = 0; i < ds->n; i++) fprintf(fp3, "%d\n", startPos[i]); */
  /* fclose(fp3); */

  /* fp3 = fopen("width.txt", "w"); */
  /* for(i = 0; i < m->mode; i++) fprintf(fp3, "%d\n", widths[i]); */
  /* fclose(fp3); */

  /////////////////////////////////////////////


  j1 = 0;

  /* Update model by maximizing posterior parameters */

  while(updateBestModel(ds, m, background, labels, startPos, j1) && j1 < 1000){
    getMotifCount(m, ds, startPos, labels);
    maxLikelihood = calculateLikelihood(m, ds, labels, startPos, background, featureCounts);
    if(fp != NULL) fprintf(fp, "%lf\n", maxLikelihood);
    /* fprintf(fp2, "%lf\n", calculateLikelihood(m, ds, lpc, spc, background, featureCounts)); */

    j1++;
  }

  /* fclose(fp2); */
  
  if(fp != NULL) fclose(fp);
  copyLabels(m->mWidth, widths, m->mode);

  double s;
  for(i = 0; i < m->mode; i++){
    if((m->t)[i] == 0) continue;
    m1 = (m->motifs)[i].motif;
    s = 0;
    for(j = 0; j < (m->mWidth)[i]; j++){
      for(k = 0; k < m->featureValues; k++){
	if(((double)(m1->modeMotifCount)[k] / (m->t)[i]) <= 0.00001) continue;
	s = s - ((double)(m1->modeMotifCount)[k] / (m->t)[i]) * log((double)(m1->modeMotifCount)[k] / (m->t)[i]);
      }
      m1 = m1->next;
    }
  }

  /* Save P(X_i | mode) for all X_i's and all modes  */
  fp = fopen(likelihoodInfoFile, "w");
  for(i = 0; i < m->mode; i++) fprintf(fp, "%lf\t", ((m->t)[i] + m->alpha) / (m->n + m->mode*m->alpha));
  fprintf(fp, "\n");
  for(i = 0; i < ds->n; i++)
    {
      for(j = 0; j < m->mode; j++)
  	fprintf(fp, "%lf\t", likelihoodXi(m, ds, j, startPos[i], i, background[i]));
      fprintf(fp, "\n");
    }
  fclose(fp);
  
				 
  /* Save values in trainOut structure and free the rest */
  for(i = 0; i < ds->n; i++) free(featureCounts[i]);
  free(featureCounts);
  for(i = 0; i < ds->n; i++) free(featureCountsCopy[i]);
  free(featureCountsCopy);

  for(i = 0; i < ds->n; i++) free(labelcounts[i]);
  free(lpc);
  free(spc);
  free(queue);
  free(x);
  free(labelcounts);
  free(modeLikes);
  freeModel(m);
  to->labels = labels;
  to->startPos = startPos;
  to->motifWidth = widths;
  to->likelihood = maxLikelihood;

  return to;

}
