#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<gsl/gsl_fit.h>
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "motifops.h"
#include "traindata.h"

FILE *fp;
int citer;

double linearSlope(queueStruct *q, double xsum, double ysum, double xxsum, double xysum, int n){
  double a1;
  a1 = (n * xysum - xsum * ysum) / (n * xxsum - xsum * xsum);
  return a1;
}

double corrCoef(int *x, double *y, int n){
  double *xx, *yy, *xy, xsum, ysum, xxsum, yysum, xysum, num, denom;
  int i;
  xx = (double*)malloc(sizeof(double)*n);
  if(!xx) printMessages(0, NULL);
  yy = (double*)malloc(sizeof(double)*n);
  if(!yy) printMessages(0, NULL);
  xy = (double*)malloc(sizeof(double)*n);
  if(!xy) printMessages(0, NULL);
  xsum = 0; ysum = 0; xxsum = 0; yysum = 0; xysum = 0;
  for(i = 0; i < n; i++){
    xy[i] = (double)x[i] * y[i];
    xx[i] = (double)x[i] * x[i];
    yy[i] = y[i] * y[i];
    xsum = xsum + (double)x[i];
    ysum = ysum + y[i];
    xxsum = xxsum + xx[i];
    yysum = yysum + yy[i];
    xysum = xysum + xy[i];
  }
  num = n*xysum - xsum * ysum;
  denom = (n * xxsum - xsum * xsum) * (n * yysum - ysum * ysum);
  free(xx);
  free(yy);
  free(xy);
  //  printf("Numerator: %lf\t Denominator: %lf\n", num, denom);
  if(denom < 0.00001) return 0;
  return num/sqrt(denom);
}

double avgCorrCoef(model *m, dataSet *ds, int *labels, int *startPos, double **background, int **labelcounts){
  double *values, corrcoef, tmp;
  int i, j;
  values = (double*)malloc(sizeof(double)*m->mode);
  if(!values) printMessages(0, NULL);
  corrcoef = 0;
  for(i = 0; i < ds->n; i++){
    addRemoveDataPoint(m, ds, labels, startPos, i, -1);
    for(j = 0; j < m->mode; j++) values[j] = likelihoodXi(m, ds, j, startPos[i], i, background[i]);
    corrcoef = corrcoef + corrCoef(labelcounts[i], values, m->mode);
    addRemoveDataPoint(m, ds, labels, startPos, i, 1);
  }
  free(values);
  return corrcoef/ds->n;
}

void enqueue(queueStruct *queue, int *front, int *rear, int *isFull, int size, double value, double *sum, double *sumxy){
  if(*isFull == size) dequeue(queue, front, rear, size, sum, sumxy);
  if(*front == -1) *front = *rear = 0;
  else *rear = (*rear + 1)%size;
  queue[*rear].y = value;
  *sum = *sum + value;
  if(sumxy != NULL) *sumxy = *sumxy + value * (*rear + 1);
  if(*isFull < size) *isFull = *isFull + 1;
  /* int i; */
  /* printf("Queue sum: %lf, Average: %lf\n", *sum, *sum/size); */
  /* printf("Queue values:\n"); */
  /* for(i = 0; i < *isFull; i++) printf("%lf\t", queue[i]); */
  /* printf("\n"); */
}

void dequeue(queueStruct *queue, int *front, int *rear, int size, double *sum, double *sumxy){
  *sum = *sum - queue[*front].y;
  if(sumxy != NULL) *sumxy = *sumxy - queue[*front].y * (*front + 1);
  if(*front == *rear) *front = *rear = -1;
  else *front = (*front + 1)%size;
}

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

int bestPosteriorParameters(dataSet *ds, model *m, double *background, int *labels, int *startPos, int index, int *t, int n, int fv){
  int i, j, k, c, id, posID, modeID;
  double s, s1, max, *priors;
  max = -1;
  id = 0;
  c = 0;
  posID = 0;
  modeID = 0;
  for(i = 0; i < m->mode; i++){
    s1 = (double)(t[i] + m->alpha)/(n + m->mode*m->alpha);
    /* j = 0; */
    /* while(j < (ds->features)[index]){ */
    /*   if((ds->lookAhead)[index][j] < (m->mWidth)[i]){ */
    /* 	s = -1; */
    /* 	j = j + (ds->lookAhead)[index][j] + 1; */
    /*   } */
    /*   else{ */
    /* 	s = likelihoodXi(m, ds, i, j, index, background); */
    /* 	if(s < 0.00001) s = -1; */
    /* 	else s = s * s1; */
    /* 	j++; */
    /*   } */
    /*   if(s > max){ */
    /* 	max = s; */
    /* 	id = c; */
    /*   } */
    /*   c++; */
    /* } */

    for(j = 0; j < (ds->features)[index]; j++){
      s = likelihoodXi(m, ds, i, j, index, background);
      if(s < 0.00001) s = -1;
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

int* getStartIndices(model *m, dataSet *ds, int mode, int index, unsigned int *seed, int flag, double *background){
  int i;
  int *indices;
  double *values;
  indices = (int*)malloc(sizeof(int)*(ds->features)[index]);
  if(!indices) printMessages(0, NULL);
  values = (double*)malloc(sizeof(double)*(ds->features)[index]);
  if(!values) printMessages(0, NULL);
  for(i = 0; i < (ds->features)[index]; i++){
    indices[i] = i;
    values[i] = likelihoodXi(m, ds, mode, i, index, background);
  }
  bubbleSort(values, indices, (ds->features)[index]);
  //  printf("Index: %d, Indices: %d, %d, %d, %d, Values: %lf, %lf, %lf, %lf\n", index, indices[0], indices[1], indices[2], indices[3], values[0], values[1], values[2], values[3]);

  free(values);
  return indices;
}

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

    /* printf("S_1: %lf\n", s); */
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
    /* printf("S_3: %lf\n", s); */
    m1 = (m->motifs)[labels[i]].motif;
    for(j = 0; j < (m->mWidth)[labels[i]]; j++){
      s = s + log((m1->modeMotifCount)[(ds->data)[i][j + startPos[i]]] + 0.5) - log((m->t)[labels[i]] + m->featureValues*0.5);
      m1 = m1->next;
    }
    for(j = 0; j < startPos[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]);
    for(j = startPos[i] + (m->mWidth)[labels[i]]; j < (ds->features)[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]);

    /* printf("S_2: %lf\n", s); */

  }
  s1 = s1 * (0.5 - 1);
  s = s + s1;
  return s;
}

/* double calculateFullPosterior(dataSet *ds, model *m, double **background, int *labels, int *startPos){ */
/*   double s, s1; */
/*   int i, j, k; */
/*   motifStruct *m1; */
/*   s = 0; */
/*   s1 = 0; */
/*   for(i = 0; i < m->mode; i++){ */
/*     s = s + ((m->t)[i])*(log((m->t)[i] + m->alpha) - log(m->n + m->mode*m->alpha)); */
/*     s = s - m->lambda*(m->mWidth)[i]; */
/*     m1 = (m->motifs)[i].motif; */
/*     for(j = 0; j < (m->mWidth)[i]; j++){ */
/*       for(k = 0; k < m->featureValues; k++){ */
/* 	s1 = s1 + log((m1->modeMotifCount)[k] + 0.5) - log((m->t)[i] + m->featureValues*0.5); */
/*       } */
/*       m1 = m1->next; */
/*     } */

/*     /\* printf("S_1: %lf\n", s); *\/ */
/*   } */
/*   for(i = 0; i < ds->n; i++){ */
/*     if(startPos[i] == -1){ */
/*       for(j = 0; j < (ds->features)[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]); */
/*       if(m->zoops > 0 && m->zoops < 1) s = s + log(m->zoops); */
/*       else s = s - log((ds->features)[i]); */
/*       continue; */
/*     } */
/*     if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log((ds->features)[i]); */
/*     else s = s - log((ds->features)[i]); */
/*     /\* printf("S_3: %lf\n", s); *\/ */
/*     m1 = (m->motifs)[labels[i]].motif; */
/*     for(j = 0; j < (m->mWidth)[labels[i]]; j++){ */
/*       s = s + log((m1->modeMotifCount)[(ds->data)[i][j + startPos[i]]] + 0.5) - log((m->t)[labels[i]] + m->featureValues*0.5); */
/*       m1 = m1->next; */
/*     } */
/*     for(j = 0; j < startPos[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]); */
/*     for(j = startPos[i] + (m->mWidth)[labels[i]]; j < (ds->features)[i]; j++) if(background[i][j] > 0.00001) s = s + log(background[i][j]); */

/*     /\* printf("S_2: %lf\n", s); *\/ */

/*   } */
/*   s1 = s1 * (0.5 - 1); */
/*   s = s + s1; */
/*   return s; */
/* } */

double calculateLikelihood(model *m, dataSet *ds, int *lp, int *startPos, double **background, int **featureCounts, int minWidth){
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
      /* if(m->zoops == 0) lxi = lxi - log(featureCounts[i][(m->mWidth)[lp[i]] - minWidth]); */
      /* else if(m->zoops == 1) lxi = lxi - log(featureCounts[i][(m->mWidth)[lp[i]] - minWidth] + 1); */
      /* else lxi = lxi + log(1 - m->zoops) - log(featureCounts[i][(m->mWidth)[lp[i]] - minWidth]); */
      s = s + lxi;
    }
    else if(startPos[i] == -1 && m->zoops != 0) s = s + log(m->zoops);
    else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][lp[i]]);
    /* else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][(m->mWidth)[lp[i]] - minWidth]); */

    /* if(lxi != 0) s = s + log(lxi); */
    /* if(m->zoops > 0 && m->zoops < 1){ */
    /*   if(startPos[i] == -1) s = s + log(m->zoops); */
    /*   else s = s + log(1 - m->zoops) - log((ds->features)[i]); */
    /* } */
    /* else s = s - log((ds->features)[i]); */
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

/* double calculateLikelihoodFaster(model *m, dataSet *ds, int *lp, int *startPos, double **background, int oldLabel, int oldStart, int index, double likelihood){ */
/*   int i; */
/*   for(i =    */
/* } */

double calculateLikelihoodMode(model *m, dataSet *ds, int *labels, int *startPos, double **background, int **featureCounts, int mode, int minWidth){
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
      /* if(m->zoops == 0) s1 = s1 - log(featureCounts[i][(m->mWidth)[mode] - minWidth]); */
      /* else if(m->zoops == 1) s1 = s1 - log(featureCounts[i][(m->mWidth)[mode] - minWidth] + 1); */
      /* else s1 = s1 + log(1 - m->zoops) - log(featureCounts[i][(m->mWidth)[mode] - minWidth]); */
      s = s + s1;
    }
    else if(startPos[i] == -1 && m->zoops != 0) s = s + log(m->zoops);
    else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][mode]);
    /* else if(m->zoops > 0 && m->zoops < 1) s = s + log(1 - m->zoops) - log(featureCounts[i][(m->mWidth)[mode] - minWidth]); */

      
    /* if(s1 != 0) s = s + log(s1); */
    /* if(m->zoops > 0 && m->zoops < 1){ */
    /*   if(startPos[i] == -1) s = s + log(m->zoops); */
    /*   else s = s + log(1 - m->zoops) - log((ds->features)[i]); */
    /* } */
    /* else s = s - log((ds->features)[i]); */
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

double likelihoodXi(model *m, dataSet *ds, int mode, int start, int index, double *background){
  int i;
  double s, s1;
  motifStruct *m1;
  s = 1;
  s1 = 1;
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

double likelihoodXiFaster(model *m, dataSet *ds, int mode, int start, int index, double *background){
  int i;
  double s, s1;
  motifStruct *m1;
  s = 1;
  s1 = 1;
  m1 = (m->motifs)[mode].motif;
  for(i = 0; i < (m->mWidth)[mode]; i++){
    s = s * ((m1->modeMotifCount)[(ds->data)[index][i + start]] + 0.5) / ((m->t)[mode] + m->featureValues*0.5);
    m1 = m1->next;
    if(background[i + start] <= 0.00001) continue;
    s = s / (background[i + start]);
  }
  return s;
}

int sampleNewLabel(model *m, dataSet *ds, int *startPos, int index, unsigned int *seed, int flag, double *background, int *featureCounts, int minWidth){
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
      /* if(m->zoops == 0) values[i] = values[i] / featureCounts[(m->mWidth)[i] - minWidth]; */
      /* else if(m->zoops == 1) values[i] = values[i] / (featureCounts[(m->mWidth)[i] - minWidth] + 1); */
      /* else values[i] = values[i] * (1 - m->zoops) / featureCounts[(m->mWidth)[i] - minWidth]; */
    }
  }
  if(flag == 1) j = arrayMax(values, m->mode);
  else{
    j = sample(values, m->mode, ((double)rand_r(seed))/(RAND_MAX));
    //    j = sample(values, m->mode, ((double)rand())/(RAND_MAX));
    //    j = sample(values, m->mode, erand48(seed));
  }
  free(values);
  return j;
}

/* int sampleStartPosFast(model *m, dataSet *ds, int mode, int index, int start, unsigned int *seed, int flag, double *background, int *indices){ */
/*   double *values, sum, v, currProb; */
/*   int i, j, count; */
/*   double *values1; */
/*   sum = 0; */
/*   count = 0; */
/*   values = (double*)malloc(sizeof(double)*((ds->features)[index])); */
/*   if(!values) printMessages(0, NULL); */
/*   values1 = (double*)malloc(sizeof(double)*((ds->features)[index])); */
/*   if(!values1) printMessages(0, NULL); */
/*   for(i = 0; i < (ds->features)[index]; i++){ */
/*     values[i] = likelihoodXi(m, ds, mode, indices[i], index, background); */
/*     if(values[i] == 0) count++; */
/*     sum = sum + values[i]; */
/*   } */

/*   copyLabels(values, values1, (ds->features)[index]); */
/*   for(i = 0; i < (ds->features)[index]; i++) values[i] = values[i]/sum; */


/*   /\* printf("Values\n"); *\/ */
/*   /\* for(i = 0; i < (ds->features)[index]; i++){ *\/ */
/*   /\*   printf("%lf\t", values[i]); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */
/*   /\* printf("Indices\n"); *\/ */
/*   /\* for(i = 0; i < (ds->features)[index]; i++){ *\/ */
/*   /\*   printf("%d\t", indices[i]); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */

/*   if(count == (ds->features)[index]) j = -1; */


/*   j = 0; */
/*   v = ((double)rand())/(RAND_MAX); */
/*   /\* printf("V: %lf, sum: %lf\n",v, sum); *\/ */
/*   //  v = v * sum; */
/*   currProb = values[0]; */
/*   /\* printf("Starting prob: %lf\n", currProb); *\/ */
/*   while(v > currProb){ */
/*     j++; */
/*     currProb = currProb + values[j]; */
/*   } */


/*   /\* j = sample(values1, (ds->features)[index], ((double)rand())/(RAND_MAX)); *\/ */

/*   /\* printf("J: %d, start: %d, v: %lf\n", j, start, v); *\/ */
/*   printf("Found in %d\n", j); */
/*   if(j != start){  */
/*     /\* bubbleSort(values, indices, (ds->features)[index]); *\/ */
/*     bubble(values, indices, j, (ds->features)[index]); */
/*     bubble(values, indices, start, (ds->features)[index]); */
/*     /\* printf("Sorted indices: %d, %d, %d, %d\n", indices[0], indices[1], indices[2], indices[3]); *\/ */
/*     /\* printf("Sorted Values: %lf, %lf, %lf, %lf\n", values[0], values[1], values[2], values[3]); *\/ */
/*   } */
/*   free(values); */
/*   free(values1); */
/*   return j; */
/* } */

int sampleStartPosFast(model *m, dataSet *ds, int mode, int index, int start, unsigned int *seed, int flag, double *background){
  double *values, sum;
  int i, j, count, *positions, size;
  size = 5;
  sum = 0;
  values = (double*)malloc(sizeof(double)*size);
  if(!values) printMessages(0, NULL);
  positions = nRandomPos(size, (ds->features)[index], start, seed);
  for(i = 0; i < size; i++){
    values[i] = likelihoodXi(m, ds, mode, positions[size], index, background);
    if(values[i] == 0) count++;
    sum = sum + values[i];
  }
  if(count == (ds->features)[index]) j = -1;
  else j = sample(values, size, ((double)rand())/(RAND_MAX));
  if(j != -1) j = positions[j];
  free(values);
  free(positions);
  /* printf("Index: %d, start: %d, Positions:\n", index, start); */
  /* for(i = 0; i < 5; i++) printf("%d\t", positions[i]); */
  /* printf("\n"); */
  return j;
}

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
  /* printf("Values:\n"); */
  /* for(i = 0; i < (ds->features)[index]; i++) */
  /*   printf("%lf\t", values[i]); */
  /* printf("\n"); */

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
    //    j = sample(values, (ds->features)[index] + (m->zoops > 0), erand48(seed));
    //    j = sample(values, (ds->features)[index] + (m->zoops > 0), ((double)rand())/(RAND_MAX));
    j = sample(values, (ds->features)[index] + (m->zoops > 0), ((double)rand_r(seed))/(RAND_MAX));
  }
  free(values);
  if(j == (ds->features)[index]) j = -1;
  return j;
}

int sampleMotifWidthLeft(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int minWidth, unsigned int *seed){
  int i, j;
  double *values;
  values = (double*)malloc(sizeof(double)*3);
  if(!values) printMessages(0, NULL);
  values[0] = 1;
  values[1] = motifRemoveLeftScore(ds, m, background, labels, startPos, mode);
  //  values[1] = 0;
  values[2] = motifAddLeftScore(ds, m, background, labels, startPos, mode);
  /* if((m->mWidth)[mode] == m->minWidth) values[1] = 0; */
  /* else values[1] = motifRemoveLeftScore(ds, m, background, labels, startPos, mode); */
  /* //  printf("Value1: %lf\n", values[1]); */
  /* if((m->mWidth)[mode] == m->maxWidth) values[2] = 0; */
  /* else values[2] = motifAddLeftScore(ds, m, background, labels, startPos, mode); */
  /* if(values[1] == values[2] && values[0] == values[1]) return 0; */
  /* if(values[1] < 0.00001 && values[2] < 0.00001) return 0; */
  if((m->mWidth)[mode] == minWidth) values[1] = 0;
  //  if((m->mWidth)[mode] == maxWidth) values[2] = 0;
  /* printf("Checking in left:\n"); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */

  //  printf("Value2: %lf\n", values[2]);
  /* printf("Old width: %d\n", (m->mWidth)[mode]); */
  /* printf("Old start: "); */
  /* i = 0; */
  /* j = 0; */
  /* while(i < 5 && j < ds->n){ */
  /*   if(labels[j++] != mode) continue; */
  /*   printf("%d\t", startPos[j - 1]); */
  /*   i++; */
  /* } */
  /* printf("\n"); */
  /* printf("Values Left:\n"); */
  /* for(i = 0; i < 3; i++) printf("%lf\t", values[i]); */
  /* printf("\n"); */
  i = sample(values, 3, ((double)rand_r(seed))/(RAND_MAX));
  //  i = sample(values, 3, ((double)rand())/(RAND_MAX));
  //  i = sample(values, 3, erand48(seed));
  if(i == 1) motifLeftDecrease(ds, m, background, labels, startPos, mode);
  else if(i == 2) motifLeftIncrease(ds, m, background, labels, startPos, mode);

  /* printf("Checking in left 2: %d\n", i); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */

  /* printf("In Left 1\n"); */
  /* printf("New Width: %d\n", (m->mWidth)[mode]); */
  /* printf("New start: "); */
  /* i = 0; */
  /* j = 0; */
  /* while(i < 5 && j < ds->n){ */
  /*   if(labels[j++] != mode) continue; */
  /*   printf("%d\t", startPos[j - 1]); */
  /*   i++; */
  /* } */
  /* printf("\n"); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */
  /* printf("In Left 2\n"); */
  free(values);
  return i;
}

int sampleMotifWidthRight(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, int minWidth, unsigned int *seed){
  int i;
  double *values;
  values = (double*)malloc(sizeof(double)*3);
  if(!values) printMessages(0, NULL);
  values[0] = 1;
  values[1] = motifRemoveRightScore(ds, m, background, labels, startPos, mode);
  values[2] = motifAddRightScore(ds, m, background, labels, startPos, mode);  
  if((m->mWidth)[mode] == minWidth) values[1] = 0;
  //  if((m->mWidth)[mode] == maxWidth) values[2] = 0;
  /* if((m->mWidth)[mode] == m->minWidth) values[1] = 0; */
  /* else values[1] = motifRemoveRightScore(ds, m, background, labels, startPos, mode); */
  /* if((m->mWidth)[mode] == m->maxWidth) values[2] = 0; */
  /* else values[2] = motifAddRightScore(ds, m, background, labels, startPos, mode);   */
  /* if(values[1] == values[2] && values[0] == values[1]) return 0; */
  /* if(values[1] < 0.00001 && values[2] < 0.00001) return 0; */

  /* printf("Checking in right:\n"); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */


  //  printf("Value1: %lf\n", values[1]);
  //  printf("Value2: %lf\n", values[2]);
  /* printf("Values Right: total: %d\n", (m->t)[mode]); */
  /* for(i = 0; i < 3; i++) printf("%lf\t", values[i]); */
  /* printf("\n"); */
  i = sample(values, 3, ((double)rand_r(seed))/(RAND_MAX));
  //  i = sample(values, 3, erand48(seed));
  //  i = sample(values, 3, ((double)rand())/(RAND_MAX));
  if(i == 2) motifRightIncrease(ds, m, background, labels, startPos, mode);
  else if(i == 1) motifRightDecrease(ds, m, background, labels, startPos, mode);

  /* printf("Checking in right 2: %d\n", i); */
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */

  free(values);
  /* if(i == 1) printf("INCREASED\n\n"); */
  /* else if(i == 2) printf("DECREASED\n\n"); */
  return i;
}

int shiftMotif(dataSet *ds, model *m, double **background, int *labels, int *startPos, int mode, unsigned int *seed){
  int i;
  double *values;
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */
  values = (double*)malloc(sizeof(double)*3);
  if(!values) printMessages(0, NULL);
  values[0] = 1;
  values[1] = motifShiftLeftScore(ds, m, background, labels, startPos, mode);
  values[2] = motifShiftRightScore(ds, m, background, labels, startPos, mode);
  i = sample(values, 2, ((double)rand())/(RAND_MAX));
  if(i == 1) motifShiftLeft(ds, m, background, labels, startPos, mode);
  else if(i = 2) motifShiftRight(ds, m, background, labels, startPos, mode);
  free(values);
  /* checkMotifCount(ds, m, labels, startPos, mode, (m->mWidth)[mode]); */
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

trainOut* trainData(dataSet *ds, int mode, float alpha, float lambda, double zoops, unsigned int seed, double **background, int *mWidth, int minWidth, char *filename){
  double maxLikelihood, tmpLikelihood, *modeLikes, tl;
  int i, j, j1, k, oldLabel, oldStart, flag, count, iterLike;
  FILE *fp1, *fp3, *fp4, *fp5;
  int *lpc, *spc;
  int **istart;
  int *oldlabels;
  motifStruct *m1;
  int *widths;
  int **labelcounts;
  double c0, c1, cov00, cov01, cov11, chisq, *x;
  int oldChangeIndex, motifWidthChange;
  int front, rear, qSize, movingCounter, movingOffset, isFull;
  double qSum, xsum, xxsum, xysum;
  queueStruct *queue;
  trainOut *to;
  int *labels, *startPos;
  int **featureCounts;
  model *m;

  featureCounts = getFeatureCounts(ds, mode, mWidth[0]);
  to = (trainOut*)malloc(sizeof(trainOut));
  if(!to) printMessages(0, NULL);
  labels = (int*)malloc(sizeof(int)*ds->n);
  if(!labels) printMessages(0, NULL);
  startPos = (int*)malloc(sizeof(int)*ds->n);
  if(!startPos) printMessages(0, NULL);

  /* featureCounts = getFeatureCounts(ds, 8, 15); */
  /* for(i = 0; i < ds->n; i++){ */
  /*   for(j = 0; j < 15 - 8 + 1; j++){ */
  /*     printf("%d\t", featureCounts[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("Length: %d\n", (ds->features)[0]); */
  /* exit(0); */


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
  


  fp = NULL;
  if(filename[0] != '0') fp = fopen(filename, "w");

  /* fp3 = fopen("EMLIKE.txt", "w"); */
  /* fp4 = fopen("movingAvg.txt", "w"); */
  /* fp5 = fopen("slopeC.txt", "w"); */

  lpc = (int*)malloc(sizeof(int)*ds->n);
  if(!lpc) printMessages(0, NULL);
  spc = (int*)malloc(sizeof(int)*ds->n);
  if(!spc) printMessages(0, NULL);
  widths = (int*)malloc(sizeof(int)*m->mode);
  if(!widths) printMessages(0, NULL);
  oldlabels = (int*)malloc(sizeof(int)*ds->n);
  if(!oldlabels) printMessages(0, NULL);

  labelcounts = (int**)malloc(sizeof(int*)*ds->n);
  if(!labelcounts) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++){
    labelcounts[i] = (int*)malloc(sizeof(int)*m->mode);
    if(!labelcounts[i]) printMessages(0, NULL);
    for(j = 0; j < m->mode; j++) labelcounts[i][j] = 0;
    labelcounts[i][labels[i]]++;
  }

  copyLabels(labels, lpc, ds->n);
  copyLabels(labels, oldlabels, ds->n);
  copyLabels(startPos, spc, ds->n);
  copyLabels(m->mWidth, widths, m->mode);

  modeLikes = (double*)malloc(sizeof(double)*m->mode);
  if(!modeLikes) printMessages(0, NULL);

  for(i = 0; i < m->mode; i++) modeLikes[i] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, i, minWidth);
  j = 0;
  i = 0;
  maxLikelihood = calculateLikelihood(m, ds, lpc, spc, background, featureCounts, minWidth);

  //  printf("Likelihood: %lf\n", maxLikelihood);
  tmpLikelihood = maxLikelihood;

  enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);

  count = 0;
  flag = 0;
  iterLike = 0;


  //  while(j < 5*ds->n){
  while(j < ds->n){
  //while(j < 10){
    //    printf("%d\n", j);
    for(i = 0; i < ds->n; i++){
      //      printf("i: %d\n", i);
      oldLabel = lpc[i];
      oldStart = spc[i];
      addRemoveDataPoint(m, ds, lpc, spc, i, -1);
      
      if(lpc[i] != -1) spc[i] = sampleStartPos(m, ds, lpc[i], i, &seed, 0, background[i]);
      if(spc[i] != -1){
      	/* if((m->t)[lpc[i]] > 20) */ lpc[i] = sampleNewLabel(m, ds, spc, i, &seed, 0, background[i], featureCounts[i], minWidth);
	labelcounts[i][lpc[i]]++;
      }

      addRemoveDataPoint(m, ds, lpc, spc, i, 1);

      if(oldLabel == lpc[i] && oldStart == spc[i]){
	/////////////////////// QUEUE OPS /////////////////////
	enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
	movingCounter = (movingCounter + 1)%(movingOffset + 1);
	if(movingCounter == movingOffset && isFull == qSize){
	  movingCounter = 0;
	  /* fprintf(fp4, "%lf\n", qSum/qSize); */
	  //	  gsl_fit_linear(x, 1, queue, 1, qSize, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
	  //	  printf("GSL: %lf, calculated: %lf\n", c1, linearSlope(queue, qSize));
	  c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	  if((c1 < 0 && c1 > -0.001) || (c1 > 0 && c1 < 0.001)){
	    flag = 1;
	  }
	  /* fprintf(fp5, "%lf\n", c1); */
	}

      ///////////////////////////////////////////////////////
	continue;
      }

      modeLikes[oldLabel] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, oldLabel, minWidth);
      if(oldLabel != lpc[i]) modeLikes[lpc[i]] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, lpc[i], minWidth);
      tmpLikelihood = 0;
      for(k = 0; k < m->mode; k++) tmpLikelihood = tmpLikelihood + modeLikes[k];
      if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
      //      fprintf(fp3, "%lf\n", EMLike(m, ds, lpc, spc, background));

      /////////////////////// QUEUE OPS /////////////////////
      
      enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
      movingCounter = (movingCounter + 1)%(movingOffset + 1);
      if(movingCounter == movingOffset && isFull == qSize){
	movingCounter = 0;
	//	fprintf(fp4, "%lf\n", qSum/qSize);
	/* gsl_fit_linear(x, 1, queue, 1, qSize, &c0, &c1, &cov00, &cov01, &cov11, &chisq); */
	c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	if((c1 < 0 && c1 > -0.001) || (c1 > 0 && c1 < 0.001)){
	  flag = 1;
	}
	//	fprintf(fp5, "%lf\n", c1);
      }

      ///////////////////////////////////////////////////////

      if(tmpLikelihood > maxLikelihood){
	copyLabels(lpc, labels, ds->n);
	copyLabels(spc, startPos, ds->n);
	copyLabels(m->mWidth, widths, m->mode);
	maxLikelihood = tmpLikelihood;
	iterLike = j;
      }
    }


    j++;
    count++;
    /* if(j > 200){ */
    /*   printf("Average: %lf\n", avgCorrCoef(m, ds, lpc, spc, background, labelcounts)); */
    /*   if(avgCorrCoef(m, ds, lpc, spc, background, labelcounts) > 0.9) break; */
    /* } */
    //    if(j > 0) continue;

    if(j%(ds->n/10) != 0) continue;
    /* //    if(flag == 0) continue; */
    /* flag = 0; */

    ////////////////////////////////////////////////////////////////////
    if(j < ds->n/10){
      oldChangeIndex = -1;
      motifWidthChange = 0;
    }
    if(oldChangeIndex > 0 && motifWidthChange >= 2 && ((c1 < 0 && c1 > -0.001) || (c1 > 0 && c1 < 0.001))) break;
    /////////////////////////////////////////////////////////////////////

    /* /\* else if(j <= 100){ *\/ */
    /* /\*   oldChangeIndex = -1; *\/ */
    /* /\*   motifWidthChange = 0; *\/ */
    /* /\* } *\/ */
    //    if(j%20 != 0) continue;
    oldChangeIndex = j;
    count = 0;
    for(i = 0; i < m->mode; i++){
      while(1){
	flag = 0;
	if((m->t)[i] < 10){ 
	  break;
	}
	else{
	  k = sampleMotifWidthLeft(ds, m, background, lpc, spc, i, minWidth, &seed);
	  if(k == 2) flag = 1;
	  j1 = sampleMotifWidthRight(ds, m, background, lpc, spc, i, minWidth, &seed);
	  if(j1 == 2) flag = 1;
	  k = k + j1;
	}
	if(k == 0){
	  motifWidthChange++;
	  break;
	}
	//	printf("Sampling for mode %d, width: %d, flag: %d\n", i, (m->mWidth)[i], flag);
	motifWidthChange = 0;
	modeFeatureCount(ds, featureCounts, (m->mWidth)[i], i);
	tmpLikelihood = calculateLikelihood(m, ds, lpc, spc, background, featureCounts, minWidth);
	
	enqueue(queue, &front, &rear, &isFull, qSize, tmpLikelihood, &qSum, &xysum);
	movingCounter = (movingCounter + 1)%(movingOffset + 1);
	if(movingCounter == movingOffset && isFull == qSize){
	  movingCounter = 0;
	  c1 = linearSlope(queue, xsum, qSum, xxsum, xysum, qSize);
	}
	modeLikes[i] = calculateLikelihoodMode(m, ds, lpc, spc, background, featureCounts, i, minWidth);
	if(fp != NULL) fprintf(fp, "%lf\n", tmpLikelihood);
	if(tmpLikelihood > maxLikelihood){
	  maxLikelihood = tmpLikelihood;
	  iterLike = j;
	  copyLabels(lpc, labels, ds->n);
	  copyLabels(spc, startPos, ds->n);
	  copyLabels(m->mWidth, widths, m->mode);
	}
	//	printf("flag: %d, motif width: %d\n", flag, (m->mWidth)[i]);
	if(flag == 0) break;
      }
    }
  }
  copyLabels(widths, m->mWidth, m->mode);
  freeMotifs(m->motifs, m->mode);
  m->motifs = initializeMotifs(m->mWidth, m->mode);
  getMotifCount(m, ds, startPos, labels);

  /* printf("Maximum likelihood: %lf\n", maxLikelihood); */
  /* printf("Iteration: %d\n", iterLike); */
  /* printf("Full posterior: %lf\n", calculateFullPosterior(ds, m, background, labels, startPos)); */

  /* printf("Average: %lf\n", avgCorrCoef(m, ds, labels, startPos, background, labelcounts)); */


  //  maxLikelihood = calculateFullPosterior(ds, m, background, labels, startPos, featureCounts);
  //  maxLikelihood = EMLike(m, ds, labels, startPos, background);
  j1 = 0;

  while(updateBestModel(ds, m, background, labels, startPos, j1) && j1 < 1000){
    getMotifCount(m, ds, startPos, labels);
    maxLikelihood = calculateLikelihood(m, ds, labels, startPos, background, featureCounts, minWidth);
    if(fp != NULL) fprintf(fp, "%lf\n", maxLikelihood);
    //    maxLikelihood = calculateFullPosterior(ds, m, background, labels, startPos, featureCounts);
    j1++;
  }

  if(fp != NULL) fclose(fp);
  /* fclose(fp3); */
  /* fclose(fp4); */
  /* fclose(fp5); */

  /* printf("Widths:\n"); */
  /* for(i = 0; i < m->mode; i++) printf("%d\t", (m->mWidth)[i]); */
  /* printf("\n"); */

  copyLabels(m->mWidth, widths, m->mode);

  /* printf("Last likelihood: %lf\n", tmpLikelihood); */

  /* printf("Size:\n"); */
  /* for(i = 0; i < m->mode; i++) */
  /*   printf("%d\t", (m->t)[i]); */
  /* printf("\n"); */

  double s;
  //  printf("Entropy:\n");
  for(i = 0; i < m->mode; i++){
    if((m->t)[i] == 0) continue;
    //    printf("Mode %d:\t", i);
    m1 = (m->motifs)[i].motif;
    s = 0;
    for(j = 0; j < (m->mWidth)[i]; j++){
      for(k = 0; k < m->featureValues; k++){
	if(((double)(m1->modeMotifCount)[k] / (m->t)[i]) <= 0.00001) continue;
	s = s - ((double)(m1->modeMotifCount)[k] / (m->t)[i]) * log((double)(m1->modeMotifCount)[k] / (m->t)[i]);
      }
      m1 = m1->next;
    }
    //    printf("%lf\n", s);
  }
  for(i = 0; i < ds->n; i++) free(featureCounts[i]);
  free(featureCounts);
  for(i = 0; i < ds->n; i++) free(labelcounts[i]);
  free(lpc);
  free(spc);
  free(queue);
  free(x);
  free(oldlabels);
  free(labelcounts);
  free(modeLikes);
  freeModel(m);
  //  printf("Maximum likelihood: %lf\n", maxLikelihood);

  to->labels = labels;
  to->startPos = startPos;
  to->motifWidth = widths;
  to->likelihood = maxLikelihood;

  return to;

}
