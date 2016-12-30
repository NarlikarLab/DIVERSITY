#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "linkedListOps.h"
#include "modelops.h"
#include "messages.h"

/* Read FASTA file into a structure */
dataSet* getData(char *s, char *outFile, int rev, int mask){
  int *features, n, tmpfeatures;
  FILE *dataFile, *fo;
  char c;
  dataSet *ds;
  int flag, ot, i, j, k, tu;
  dataFile = fopen(s, "r");
  if(!dataFile) printMessages(4, s);
  n = 0;
  tmpfeatures = 0;
  flag = 0;
  /* Check if it is a valid Fasta file */
  while((c = fgetc(dataFile)) != EOF)
    if(c == '>' && flag == 0){ 
      n++;
      flag = 1;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
  rewind(dataFile);
  if(n == 0) printMessages(1, s);
  features = (int*)malloc(sizeof(int)*n);
  if(!features) printMessages(0, NULL);
  for(i = 0; i < n; i++) features[i] = 0;
  ot = 0;
  c = fgetc(dataFile);
  if(c != '>') printMessages(1, s);
  rewind(dataFile);
  tu = -1;
  flag = 0;
  i = 0;
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      if(ot == 1 && tmpfeatures == 0) printMessages(2, s);
      if(ot == 0 && tmpfeatures > 0) ot = 1;
      if(tmpfeatures > 0){
	features[i] = tmpfeatures;
	i++;
      }
      flag = 1;
      tmpfeatures = 0;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
    else if(c == '\n' && flag == 0)
      continue;
    else if(((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) && flag == 0){
      switch(c){
      case 'A': tmpfeatures++; break;
      case 'a': tmpfeatures++; break;
      case 'C': tmpfeatures++; break;
      case 'c': tmpfeatures++; break;
      case 'G': tmpfeatures++; break;
      case 'g': tmpfeatures++; break;
      case 'T': tmpfeatures++; break;
      case 't': tmpfeatures++; break;
      case 'N': tmpfeatures++; break;
      default: printf("ERROR: Invalid symbol %c in file %s\n", c, s); exit(1);
      }
    }
    if((c == 'T' || c == 't') && tu == -1 && flag == 0) tu = 0;
    else if((c == 'U' || c == 'u') && tu == -1 && flag == 0) tu = 1;
    if((c == 'T' || c == 't') && tu == 1 && flag == 0) printMessages(3, s);
    else if((c == 'U' || c == 'u') && tu == 0 && flag == 0) printMessages(3, s);
  }

  if(tmpfeatures == 0) printMessages(2, s);
  else features[i] = tmpfeatures;

  if(rev == 1)
    for(i = 0; i < n; i++)
      features[i] = 2 * features[i] + 1;
  rewind(dataFile);
  ds = (dataSet*)malloc(sizeof(dataSet));
  if(!(ds)) printMessages(0, NULL);
  ds->data = (int**)malloc(sizeof(int*)*n);
  if(!(ds->data)) printMessages(0, NULL);
  for(i = 0; i < n; i++){
    (ds->data)[i] = (int*)malloc(sizeof(int)*features[i]);
    if(!((ds->data)[i])) printMessages(0, NULL);
  }
  i = -1;
  flag = 0;
  j = 0;
  /* Save data in structure */
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      i++;
      j = 0;
      flag = 1;
    }
    else if(c == '\n' && flag == 1)
      flag = 0;
    else if(c == '\n' && flag == 0)
      continue;
    else if(((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) && flag == 0){
      switch(c){
      case 'A': (ds->data)[i][j] = 0; j++; break;
      case 'a': (ds->data)[i][j] = 0; j++; break;
      case 'C': (ds->data)[i][j] = 1; j++; break;
      case 'c': (ds->data)[i][j] = 1; j++; break;
      case 'G': (ds->data)[i][j] = 2; j++; break;
      case 'g': (ds->data)[i][j] = 2; j++; break;
      case 'T': (ds->data)[i][j] = 3; j++; break;
      case 't': (ds->data)[i][j] = 3; j++; break;
      case 'N': (ds->data)[i][j] = 4; j++; break;
      default: printf("ERROR: Invalid symbol %c in file %s\n", c, s); exit(1);
      }
      if(c > 96 && mask == 1) (ds->data)[i][j - 1] = 4;
    }
  }
  rewind(dataFile);
  fo = fopen(outFile, "w");
  flag = 0;
  ot = 0;
  while((c = fgetc(dataFile)) != EOF){
    if(c == '>' && flag == 0){
      if(ot != 0) fprintf(fo, "\n");
      else ot = 1;
      flag = 1;
    }
    else if(flag == 1 && c == '\n'){
      fprintf(fo, "\t");
      flag = 0;
    }
    else if(flag == 1){
      if(!(c <= 31 || c >=127))
	fprintf(fo, "%c", c);
    }
    else if(flag == 0 && (c == '\n' || c == ' ' || c == '\t')) continue;
    else{
      if(!(c <= 31 || c >=127)){
	if(c > 96 && mask == 1) fprintf(fo, "N");
	else fprintf(fo, "%c", c);
      }
    }
  }
  fclose(dataFile);
  fclose(fo);

  ds->n = n;
  ds->features = features;
  ds->featureValues = 4;
  ds->tu = tu;

  if(rev == 1){
    for(i = 0; i < ds->n; i++){
      (ds->data)[i][(ds->features)[i] / 2] = 4;
      j = 0;
      k = (ds->features)[i] - 1;
      while(k > ((ds->features)[i] / 2)){
	(ds->data)[i][k] = 3 - (ds->data)[i][j];
	if((ds->data)[i][k] < 0) (ds->data)[i][k] = 4;
	j++; k--;
      }
    }
  }
  ds->lookAhead = lookForN(ds);
  return ds;
}

/* Compute background  */
double** getBackground(dataSet *ds, int order){
  double **back, **counts, lsum;
  int i, j, k, l, ind, start;
  back = (double**)malloc(sizeof(double*)*ds->n);
  if(!back) printMessages(0, NULL);
  counts = (double**)malloc(sizeof(double*)*(order + 1));
  if(!counts) printMessages(0, NULL);
  for(i = 0; i <= order; i++){
    counts[i] = (double*)malloc(sizeof(double)*pow(4, i + 1));
    if(!(counts[i])) printMessages(0, NULL);
  }
  for(i = 0; i < ds->n; i++){
    back[i] = (double*)malloc(sizeof(double)*(ds->features)[i]);
    if(!(back[i])) printMessages(0, NULL);
  }
  for(i = 0; i < ds->n; i++){
    for(j = 0; j <= order; j++)
      for(k = 0; k < pow(4, j + 1); k++)
	counts[j][k] = 0;
    for(j = 0; j < (ds->features)[i]; j++){
      if((ds->data)[i][j] == 4) continue;
      for(k = 0; k <= order; k++){
	ind = getIndex((ds->data)[i], j, (ds->features)[i], k);
	if(ind == -1) continue;
	counts[k][ind]++;
      }
    }
    for(j = 0; j <= order; j++){
      for(k = 0; k < pow(4, j); k++){
	lsum = 0;
	for(l = 0; l < 4; l++)
	  lsum = lsum + counts[j][k*4 + l];
	if(lsum == 0) lsum = 1;
	for(l = 0; l < 4; l++)
	  counts[j][k*4 + l] = (double)counts[j][k*4 + l] / lsum;
      }
    }

    for(j = 0; j < (ds->features)[i]; j++){
      if((ds->data)[i][j] == 4){
	back[i][j] = -1;
	continue;
      }
      start = maxWithoutN((ds->data)[i], j, order);
      ind = getIndex((ds->data)[i], start, (ds->features)[i], j - start);
      if(ind == -1) back[i][j] = -1;
      else back[i][j] = counts[j - start][ind];
    }
  }

  for(i = 0; i <= order; i++)
    free(counts[i]);
  free(counts);
  return back;
}

/* Compute maximum. Return 0 if start position is negative */
int maxWithoutN(int *data, int pos, int order){
  int i, start;
  start = pos - order;
  if(start < 0) start = 0;
  for(i = pos; i >= start; i--)
    if(data[i] == 4) return i + 1;
  if(start < 0) return -1;
  return start;
}

/* Compute index in the given sequence array */
int getIndex(int *sequence, int start, int end, int order){
  int i, index;
  index = 0;
  for(i = 0; (i <= order && (start + i) < end); i++){
    if(!(sequence[start + i] >= 0 && sequence[start + i] < 4)) return -1;
    index = index + pow(4, order - i)*sequence[start + i];
  }
  return index;
}

/* Free structure dataSet */
void freeData(dataSet *ds){
  int i;
  for(i = 0; i < ds->n; i++){
    free((ds->data)[i]); (ds->data)[i] = NULL;
  }
  free(ds->lookAhead); ds->lookAhead = NULL;
  free(ds->features); ds->features = NULL;
  free(ds->data); ds->data = NULL;
  free(ds); ds = NULL;
}

/* Create model structure */
model* createModel(int mode, dataSet *ds, int *labels, int *startPos, float alpha, float lambda, double zoops, int *mWidth){
  int i;
  model *m = initializeModel(mode, ds->featureValues, mWidth);
  m->n = ds->n;
  m->mode = mode;
  m->featureValues = ds->featureValues;
  m->alpha = alpha;
  m->zoops = zoops;
  m->lambda = lambda;
  for(i = 0; i < mode; i++) (m->mWidth)[i] = mWidth[i];
  getMotifCount(m, ds, startPos, labels);
  return m;
}

/* Initialize model structure */
model *initializeModel(int mode, int featureValues, int *mWidth){
  model *m;
  m = (model*)malloc(sizeof(model));
  if(!m) printMessages(0, NULL);
  m->motifs = initializeMotifs(mWidth, mode);
  m->t = (int*)malloc(sizeof(int)*mode);
  if(!m->t) printMessages(0, NULL);
  m->mWidth = (int*)malloc(sizeof(int)*mode);
  if(!m->mWidth) printMessages(0, NULL);
  return m;
}

/* Free model structure */
void freeModel(model *m){
  freeMotifs(m->motifs, m->mode);
  free(m->t);
  m->t = NULL;
  free(m->mWidth);
  m->mWidth = NULL;
  free(m);
}

/* Compute motif counts based on the start position of the motifs and the motif widths */
void getMotifCount(model *m, dataSet *ds, int *startPos, int *labels){
  int i, j, k, mode;
  motifStruct *m1;
  for(i = 0; i < m->mode; i++){
    m1 = (m->motifs)[i].motif;
    for(j = 0; j < (m->mWidth)[i]; j++){
      for(k = 0; k < ds->featureValues; k++)
	(m1->modeMotifCount)[k] = 0;
      m1 = m1->next;
    }
    (m->t)[i] = 0;
  }
  m->n = 0;
  for(i = 0; i < ds->n; i++){
    mode = labels[i];
    if(startPos[i] < 0) continue;
    (m->t)[mode]++;
    m->n++;
    
    m1 = (m->motifs)[mode].motif;
    for(j = 0; j < (m->mWidth)[mode]; j++){
      (m1->modeMotifCount)[(ds->data)[i][j + startPos[i]]]++;
      m1 = m1->next;
    }
  }
}

/* Sample values for labels array */
int* sampleArrayVals(int n, int mode, int addVal, unsigned int *seed){
  int i, *lp;
  lp = (int*)malloc(sizeof(int)*n);
  if(!lp) printMessages(0, NULL);
  for(i = 0; i < n; i++) lp[i] = (rand()%mode)+addVal;
  return lp;
}

/* Sample start positions */
int* sampleStart(int n, int vals, int zoops, int notPresent, unsigned int *seed){
  int i, *lp;
  lp = (int*)malloc(sizeof(int)*n);
  if(!lp) printMessages(0, NULL);
  for(i = 0; i < n; i++){
    lp[i] = (rand()%(vals + zoops));
    if(lp[i] == vals) lp[i] = notPresent;
  }
  return lp;
}

/* Compute exp from given log values */
void logToExp(double *p, int n){
  double max;
  int i;
  max = p[0];
  for(i = 0; i < n; i++)
    max = p[i] > max ? p[i] : max;
  for(i = 0; i < n; i++)
    p[i] = pow(exp(1), p[i] - max);
}

/* Sample a value based on given ditribution */
int sample(double *p, int n, double r){
  int low, high, mid, i;
  double sum;
  low = 0;
  high = n;
  mid = n/2;
  sum = 0;
  for(i = 0; i < n; i++)
    sum = sum + p[i];
  for(i = 0; i < n; i++)
    p[i] = p[i]/sum;
  for(i = 1; i < n; i++)
    p[i] = p[i] + p[i-1];
  while(low <= high){
    mid = (low + high)/2;
    if(r == p[mid]){
      if(mid >= n) return n-1;
      return mid;
    }
    else if(r > p[mid])
      low = mid + 1;
    else if(r < p[mid])
      high = mid - 1;
  }
  if(low >= n) return n-1;
  return low;
}

/* Bubble sort */
void bubbleSort(double *values, int *indices, int n){
  int i, j, temp;
  double tempV;
  for(i = 0; i < n; i++){
    for(j = 0; j < n - 1; j++){
      if(values[j + 1] > values[j]){
	tempV = values[j];
	temp = indices[j];
	values[j] = values[j + 1];
	indices[j] = indices[j + 1];
	values[j + 1] = tempV;
	indices[j + 1] = temp;
      }
    }
  }
}


/* Return index of maximum value in array */
int arrayMax(double *p, int n){
  int i, index;
  double val;
  
  val = p[0];
  index = 0;
  for(i = 1; i < n; i++)
    if(p[i] > val){
      val = p[i];
      index = i;
    }
  return index;
}

/* Free structure trainOut */
void freeTo(trainOut *to){
  free(to->labels); to->labels = NULL;
  free(to->startPos); to->startPos = NULL;
  free(to->motifWidth); to->motifWidth = NULL;
  free(to);
}

/* Copy model to another model structure */
void copyModel(model *m, model *mc){
  int i, j, k;
  motifStruct *m1, *m2;
  mc->mode = m->mode;
  mc->featureValues = m->featureValues;
  mc->n = m->n;
  mc->alpha = m->alpha;
  mc->zoops = m->zoops;
  for(i = 0; i < mc->mode; i++){
    (mc->mWidth)[i] = (m->mWidth)[i];
    m1 = (mc->motifs)[i].motif;
    m2 = (m->motifs)[i].motif;
    for(j = 0; j < (mc->mWidth)[i]; j++){
      for(k = 0; k < mc->featureValues; k++){
	(m1->modeMotifCount)[k] = (m2->modeMotifCount)[k];
      }
      m1 = m1->next;
      m2 = m2->next;
    }
    (mc->t)[i] = (m->t)[i];
  }
}

/* Copy values from one array to another */
void copyLabels(int *lp1, int *lp2, int n){
  int i;
  for(i = 0; i < n; i++)
    lp2[i] = lp1[i];
}

/* Check if motif contains N instead of A, C, G, T */
int motifWithN(int *data, int start, int width){
  int i;
  for(i = start; i < start + width; i++)
    if(data[i] > 3) return 1;
  return 0;
}

/* Initialize labels and stat positions for a given model */
void initializeLabelStartPos(dataSet *ds, int *labels, int *startPos, int mode, int *pos, unsigned int *seed){
  int i, count;
  for(i = 0; i < ds->n; i++){
    labels[i] = rand_r(seed)%(mode);
    if((ds->features)[i] - pos[labels[i]] + 1 <= 0){
      startPos[i] = -1;
      continue;
    }
    count = 0;
    while(count < 20){
      startPos[i] = rand_r(seed)%((ds->features)[i] - pos[labels[i]] + 1);
      if(motifWithN((ds->data)[i], startPos[i], pos[labels[i]])) count++;
      else break;
    }
    if(count == 20) startPos[i] = -1;
  }
}

/* Sample n random values from given set of values */
int* nRandomPos(int posCount, int size, int lastPos, unsigned int *seed){
  int *arr, i, n, num, *numArr;
  arr = (int*)malloc(sizeof(int)*posCount);
  if(!arr) printMessages(0, NULL);
  numArr = (int*)malloc(sizeof(int)*size);
  if(!numArr) printMessages(0, NULL);
  n = size;
  for(i = 0; i < size; i++) numArr[i] = i;
  i = 0;
  while(i < (posCount - 1)){
    num = (int)((double)rand())%n;
    if(numArr[num] == lastPos) continue;
    arr[i] = numArr[num];
    numArr[num] = numArr[n - 1];
    n--;
    i++;
  }
  free(numArr);
  arr[posCount - 1] = lastPos;
  return arr;
}

/* Normalize values in array */
double* normalize(int *a, int n){
  double sum;
  double *values;
  int i;
  values = (double*)malloc(sizeof(double)*n);
  if(!values) printMessages(0, NULL);
  sum = 0;
  for(i = 0; i < n; i++){
    sum = sum + a[i];
    values[i] = (double)a[i];
  }
  for(i = 0; i < n; i++) values[i] = values[i]/sum;
  return values;
}

/* Find minimum of two values */
int min(int a, int b){
  if(a < b) return a;
  else return b;
}

int** getFeatureCounts(dataSet *ds, int mode, int initialWidth){
  int **counts, i, j, size;
  size = mode;
  counts = (int**)malloc(sizeof(int*)*ds->n);
  if(!counts) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++){
    counts[i] = (int*)malloc(sizeof(int)*size);
    if(!counts[i]) printMessages(0, NULL);
  }
  modeFeatureCount(ds, counts, initialWidth, 0);
  for(i = 0; i < ds->n; i++)
    for(j = 1; j < size; j++)
      counts[i][j] = counts[i][0];
  return counts;
}

/* Get counts for all modes */
void modeFeatureCount(dataSet *ds, int **counts, int motifWidth, int index){
  int i, j, k;
  for(i = 0; i < ds->n; i++){
    counts[i][index] = 0;
    for(j = 0; j < (ds->features)[i]; j++){
      counts[i][index]++;
      for(k = 0; k < motifWidth; k++){
	if(j + k > (ds->features)[i]){
	  counts[i][index]--;
	  break;
	}
	if((ds->data)[i][j + k] > 3){
	  counts[i][index]--;
	  break;
	}
      }
    }
  }
}

/* Check if there is N istead of A, C, G, T */
int** lookForN(dataSet *ds){
  int i, j, **lookahead, count;
  lookahead = (int**)malloc(sizeof(int*)*ds->n);
  if(!lookahead) printMessages(0, NULL);
  for(i = 0; i < ds->n; i++){
    lookahead[i] = (int*)malloc(sizeof(int)*(ds->features)[i]);
    if(!lookahead[i]) printMessages(0, NULL);
    for(j = 0; j < (ds->features)[i]; j++)
      if((ds->data)[i][j] > 3) lookahead[i][j] = 0;
    count = 0;
    for(j = (ds->features)[i] - 1; j >= 0; j--){
      if((ds->data)[i][j] < 4){
	count++;
	lookahead[i][j] = count;
      }
      else count = 0;
    }
  }
  return lookahead;
}
