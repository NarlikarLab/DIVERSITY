#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "crossValidation.h"

int* posList(int n){
  int *pos;
  unsigned int seed;
  seed = 1;
  pos = (int*)malloc(sizeof(int)*n);
  if(!pos) printMessages(0, NULL);
  nRandom1(pos, n, n, &seed);
  return pos;
}

dataSet* getTrainSubset(dataSet *ds, int subset, int kfold, int *pos){
  int startPos, endPos, i, j, tr;
  dataSet *trainSet;
  startPos = subset*ds->n/kfold;
  endPos = subset == (kfold-1) ? ds->n : (subset+1)*ds->n/kfold;
  trainSet = (dataSet*)malloc(sizeof(dataSet));
  if(!trainSet) printMessages(0, NULL);
  trainSet->n = ds->n - endPos + startPos;
  trainSet->tu = ds->tu;
  trainSet->data = (int**)malloc(sizeof(int*)*trainSet->n);
  if(!(trainSet->data)) printMessages(0, NULL);
  trainSet->features = (int*)malloc(sizeof(int)*trainSet->n);
  if(!trainSet->features) printMessages(0, NULL);
  trainSet->featureValues = ds->featureValues;
  tr = 0;
  for(i = 0; i < startPos; i++){
    (trainSet->features)[tr] = (ds->features)[pos[i]];
    (trainSet->data)[tr] = (int*)malloc(sizeof(int)*(trainSet->features)[tr]);
    if(!(trainSet->data)[tr]) printMessages(0, NULL);
    for(j = 0; j < (trainSet->features)[tr]; j++)
      (trainSet->data)[tr][j] = (ds->data)[pos[i]][j];
    tr++;
  }
  for(i = endPos; i < ds->n; i++){
    (trainSet->features)[tr] = (ds->features)[pos[i]];
    (trainSet->data)[tr] = (int*)malloc(sizeof(int)*(trainSet->features)[tr]);
    if(!(trainSet->data)[tr]) printMessages(0, NULL);
    for(j = 0; j < (trainSet->features)[tr]; j++)
      (trainSet->data)[tr][j] = (ds->data)[pos[i]][j];
    tr++;
  }
  return trainSet;
}

/* Get subset for training of given fold */
dataSet* getTestSubset(dataSet *ds, int subset, int kfold, int *pos){
  int startPos, endPos, i, j, te;
  dataSet *testSet;
  startPos = subset*ds->n/kfold;
  endPos = subset == (kfold-1) ? ds->n : (subset+1)*ds->n/kfold;
  testSet = (dataSet*)malloc(sizeof(dataSet));
  if(!testSet) printMessages(0, NULL);
  testSet->n = endPos - startPos;
  testSet->tu = ds->tu;
  testSet->data = (int**)malloc(sizeof(int*)*testSet->n);
  if(!(testSet->data)) printMessages(0, NULL);
  testSet->features = (int*)malloc(sizeof(int)*testSet->n);
  if(!testSet->features) printMessages(0, NULL);
  testSet->featureValues = ds->featureValues;
  te = 0;
  for(i = startPos; i < endPos; i++){
    (testSet->features)[te] = (ds->features)[pos[i]];
    (testSet->data)[te] = (int*)malloc(sizeof(int)*(testSet->features)[te]);
    if(!(testSet->data)[te]) printMessages(0, NULL);
    for(j = 0; j < (testSet->features)[te]; j++){
      (testSet->data)[te][j] = (ds->data)[pos[i]][j];
    }
     te++;
  }
  return testSet;
}

void nRandom1(int *pos, int posCount, int size, unsigned int *seed){
  int i, n, num, *numArr;
  numArr = (int*)malloc(sizeof(int)*size);
  if(!numArr) printMessages(0, NULL);
  n = size;
  for(i = 0; i < size; i++)
    numArr[i] = i;
  for(i = 0; i < posCount; i++){
    num = rand_r(seed)%n;
    pos[i] = numArr[num];
    numArr[num] = numArr[n-1];
    n--;
  }
  free(numArr);
}

double chooseBestStart(dataSet *testSet, model *m, double **background, int index, int mode){
  int i, j, size;
  double *p, max;
  motifStruct *ms;
  max = 0;
  size = (testSet->features)[index] - (m->mWidth)[mode];
  p = (double*)malloc(sizeof(double)*size);
  for(i = 0; i < size; i++) p[i] = 0;
  for(i = 0; i < size; i++){
    ms = (m->motifs)[mode].motif;
    for(j = 0; j < (m->mWidth)[mode]; j++){
      if((testSet->data)[index][i + j] > 3){
	p[i] = 0;
	break;
      }
      p[i] = p[i] + log((ms->modeMotifCount)[(testSet->data)[index][i + j]] + m->alpha) - log((m->t)[mode] + m->featureValues*m->alpha);
      if(background[index][i + j] > 0.00001) p[i] = p[i] - log(background[index][i + j]);
      ms = ms->next;
    }
    p[i] = p[i] + log((m->t)[mode] + m->featureValues*m->alpha) - log(m->n + m->mode*m->featureValues*m->alpha);
    for(j = 0; j < (testSet->features)[index]; j++)
      if(background[index][i + j] > 0.00001) p[i] = p[i] + log(background[index][i + j]);
    if(p[i] == 0) continue;
    if(max == 0 || p[i] > max) max = p[i];
  }
  free(p);
  return max;
}
    

//////////////////////////////////////////////// One Pos for Each Mode ///////////////////////////////////////////////////////////////
/* double cvLikelihood(dataSet *testSet, model *m, double **background){ */
/*   int i, j, maxPos; */
/*   double s, sum, max, *p; */
/*   sum = 0; */
/*   p = (double*)malloc(sizeof(double)*m->mode); */
/*   for(i = 0; i < testSet->n; i++){ */
/*     for(j = 0; j < m->mode; j++){ */
/*       p[j] = chooseBestStart(testSet, m, background, i, j); */
/*     } */
/*     max = p[0]; */
/*     maxPos = 0; */
/*     for(j = 0; j < m->mode; j++){ */
/*       if(p[j] == 0) continue; */
/*       if(max == 0 || p[j] > max){ */
/* 	max = p[j]; */
/* 	maxPos = j; */
/*       } */
/*     } */
/*     s = 0; */
/*     for(j = 0; j < m->mode; j++) */
/*       if(p[j] != 0) s = s + pow(exp(1), p[j] - max); */
/*     s = log(s) + max; */
/*     sum = sum + s; */
/*   } */
/*   free(p); */
/*   return sum; */
/* } */


////////////////////////////////////////////// One Pos /////////////////////////////////////////////////////////
/* double cvLikelihood(dataSet *testSet, model *m, double **background){ */
/*   int i, j, k, l, pos, maxFeatures, maxPos, maxM; */
/*   double s, sum, max, *p; */
/*   motifStruct *ms; */
/*   sum = 0; */
/*   maxFeatures = 0; */
/*   for(i = 0; i < testSet->n; i++) */
/*     if((testSet->features)[i] > maxFeatures) maxFeatures = (testSet->features)[i]; */
/*   p = (double*)malloc(sizeof(double)*maxFeatures*m->mode); */
/*   for(i = 0; i < testSet->n; i++){ */
/*     max = 0; */
/*     for(j = 0; j < m->mode; j++){ */
/*       for(k = 0; k < (testSet->features)[i]; k++){ */
/* 	pos = j*(testSet->features)[i] + k; */
/* 	p[pos] = 0; */
/* 	if(pos > (testSet->features)[i] - (m->mWidth)[j]) break; */
/* 	ms = (m->motifs)[j].motif; */
/* 	/\* p[pos] = log((m->t)[j] + m->featureValues*m->alpha) - log(m->n + m->mode*m->featureValues*m->alpha); // Need to add probability of feature as start *\/ */
/* 	for(l = 0; l < (m->mWidth)[j]; l++){ */
/* 	  if((testSet->data)[i][k + l] > 3){ */
/* 	    p[pos] = 0; */
/* 	    break; */
/* 	  } */
/* 	  p[pos] = p[pos] + log((ms->modeMotifCount)[(testSet->data)[i][k + l]] + m->alpha) - log((m->t)[j] + m->featureValues*m->alpha); */
/* 	  if(background[i][k + l] > 0.00001) p[pos] = p[pos] - log(background[i][k + l]); */
/* 	  ms = ms->next; */
/* 	} */
/* 	if(p[pos] == 0) continue; */
/* 	if(max == 0 || max > p[pos]) max = p[pos]; */
/* 	//	printf("P[%d]: %lf\n", pos, p[pos]); */
/*       } */
/*     } */
/*     /\* max = p[0]; *\/ */
/*     /\* maxPos = j; *\/ */
/*     /\* for(j = 1; j < m->mode*(testSet->features)[i]; j++){ *\/ */
/*     /\*   if(p[j] == 0) continue; *\/ */
/*     /\*   if(max == 0 || p[j] > max){ *\/ */
/*     /\* 	max = p[j]; *\/ */
/*     /\* 	maxPos = j; *\/ */
/*     /\*   } *\/ */
/*     /\* } *\/ */
/*     /\* s = 0; *\/ */
/*     /\* for(j = 0; j < m->mode*(testSet->features)[i]; j++) *\/ */
/*     /\*   if(p[j] != 0) s = s + pow(exp(1), p[j] - max); *\/ */
/*     /\* s = log(s) + max; *\/ */
/*     /\* sum = sum + s; *\/ */
/*     sum = sum + max; */
/*     //    printf("Seq: %d, Max: %lf, Sum: %lf\n", i, max, sum); */
/*   } */
/*   free(p); */
/*   return sum; */
/* } */


//////////////////////////////////////////////////// ALL POS PRODUCT /////////////////////////////////////////////////

/* double cvLikelihood(dataSet *testSet, model *m, double **background){ */
/*   int i, j, k, l; */
/*   double s, sum, max, *p; */
/*   motifStruct *ms; */
/*   sum = 0; */
/*   p = (double*)malloc(sizeof(double)*m->mode); */
/*   for(i = 0; i < testSet->n; i++){ */
/*     for(j = 0; j < (testSet->features)[i]; j++){ */
/*       for(k = 0; k < m->mode; k++){ */
/* 	p[k] = log((m->t)[k] + m->featureValues*m->alpha) - log(m->n + m->mode*m->featureValues*m->alpha); */
/* 	ms = (m->motifs)[k].motif; */
/* 	for(l = 0; l < (m->mWidth)[k]; l++){ */
/* 	  if((j + l) > (testSet->features)[i] - (m->mWidth)[k]){ */
/* 	    p[k] = 0; */
/* 	    break; */
/* 	  } */
/* 	  if((testSet->data)[i][j + l] > 3){ */
/* 	    p[k] = 0; */
/* 	    break; */
/* 	  } */
/* 	  p[k] = p[k] + log((ms->modeMotifCount)[(testSet->data)[i][j + l]] + m->alpha) - log((m->t)[k] + m->featureValues*m->alpha); */
/* 	  if(background[i][j + l] > 0.00001) p[k] = p[k] - log(background[i][j + l]); */
/* 	  ms = ms->next; */
/* 	} */
/*       } */
/*       max = p[0]; */
/*       for(k = 1; k < m->mode; k++){ */
/* 	if(p[k] == 0) continue; */
/* 	if(max == 0 || p[k] > max) max = p[k]; */
/*       } */
/*       s = 0; */
/*       for(k = 0; k < m->mode; k++){ */
/* 	if(p[k] != 0) s = s + pow(exp(1), p[k] - max); */
/*       } */
/*       //      printf("S: %lf\n", s); */
/*       /\* if(s > 0.00001){ *\/ */
/*       /\* 	s = log(s) + max; *\/ */
/*       /\* 	sum = sum + s; *\/ */
/*       /\* } *\/ */
/*       //      printf("s0: %lf\n", s); */
/*       if(s != 0) s = log(s) + max; */
/*       sum = sum + s; */
		   
/*       //      printf("s: %lf\n", s); */

/*     } */
/*   } */
/*   free(p); */
/*   //  printf("Returning: %lf\n", sum); */
/*   return sum; */
/* } */


///////////////////////////////////////////////////// ALL POS SUM //////////////////////////////////////////////////

double cvLikelihood(dataSet *testSet, model *m, double **background){

  int i, j, k, l, pos, *maxFeatures, modeCount;
  double s1, s2, sum, max, *p;
  motifStruct *ms;
  sum = 0;
  maxFeatures = (int*)malloc(sizeof(int)*m->mode);
  for(i = 0; i < testSet->n; i++){

    pos = 0;
    p = (double*)malloc(sizeof(double)*(testSet->features)[i]*m->mode);

    ////////////////////////////////////// Computing maxFeatures /////////////////////////////////////////
    for(j = 0; j < m->mode; j++){
      maxFeatures[j] = (testSet->features)[i];
      for(k = 0; k < (testSet->features)[i]; k++){
	for(l = 0; l < (m->mWidth)[j]; l++){
	  if(k + l > (testSet->features)[i] - (m->mWidth)[j]){
	    maxFeatures[j]--;
	    break;
	  }
	  if((testSet->data)[i][k + l] > 3){
	    maxFeatures[j]--;
	    break;
	  }
	}
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    s1 = 0;
    for(j = 0; j < (testSet->features)[i]; j++){

      modeCount = m->mode;
      for(k = 0; k < m->mode; k++){
	for(l = 0; l < (m->mWidth)[k]; l++){
	  if(j + l > (testSet->features)[i] - (m->mWidth)[k]){
	    modeCount--;
	    break;
	  }
	  if((testSet->data)[i][j + l] > 3){
	    modeCount--;
	    break;
	  }
	}
      }
      
      for(k = 0; k < m->mode; k++){

	//	s2 = ((m->t)[k] + m->featureValues*m->alpha) / (m->n + m->mode*m->featureValues*m->alpha);
	s2 = log((m->t)[k] + m->featureValues*m->alpha) - log(m->n + m->mode*m->featureValues*m->alpha);
	//	s2 = 1;
	ms = (m->motifs)[k].motif;
	for(l = 0; l < (m->mWidth)[k]; l++){
	  if(j + l > (testSet->features)[i] - (m->mWidth)[k]){
	    s2 = 0;
	    break;
	  }
	  if((testSet->data)[i][j + l] > 3){
	    s2 = 0;
	    break;
	  }
	  //	  s2 = s2 * ((ms->modeMotifCount)[(testSet->data)[i][j + l]] + m->alpha) / ((m->t)[k] + m->featureValues*m->alpha);
	  s2 = s2 + log((ms->modeMotifCount)[(testSet->data)[i][j + l]] + m->alpha) - log((m->t)[k] + m->featureValues*m->alpha);
	  //	  if(background[i][j + l] > 0.00001) s2 = s2 / background[i][j + l];
	  if(background[i][j + l] > 0.00001) s2 = s2 - log(background[i][j + l]);
	  ms = ms->next;
	}
	if(modeCount == 0 && s2 != 0) printf("ERROR!!!\n");
	if(maxFeatures[k] == 0 && s2 != 0) printf("ERROR1!!!\n");
	//	if(modeCount != 0) s2 = s2 / modeCount;
	if(modeCount != 0) s2 = s2 - log((double)modeCount);
	if(maxFeatures[k] != 0) s2 = s2 - log((double)maxFeatures[k]);
	//	if(maxFeatures[k] != 0) s2 = s2 / maxFeatures[k];
	//	s1 = s1 + s2;
	p[pos] = s2;
	pos++;
      }
    }

    max = p[0];
    for(j = 1; j < (testSet->features)[i]*m->mode; j++){
      if(p[j] == 0) continue;
      if(max == 0 || p[j] > max) max = p[j];
    }
    s1 = 0;

    for(j = 0; j < (testSet->features)[i]*m->mode; j++)
      if(p[j] != 0) s1 = s1 + pow(exp(1), p[j] - max);
    
    /* printf("S1: %lf\n", s1); */
    if(s1 > 0.00001) s1 = log(s1) + max;

    //    for(j = 0; j < (testSet->features)[i]*m->mode; j++) s1 = s1 + p[j];
    //    printf("S1: %lf\n", s1);
    sum = sum + s1;
    /* printf("Sum: %lf\n", sum); */
    //    printf("Sum: %lf\n", sum);

    free(p);
  }
  free(maxFeatures);
  return sum;
}

///////////////////////////////////////////////////////////////// ONE POS //////////////////////////////////////////////////////////////////////

/* double cvLikelihood(dataSet *testSet, model *m, double **background){ */
/*   int i, j, k, l, pos, *maxFeatures, modeCount; */
/*   double s1, s2, sum, max, *p; */
/*   motifStruct *ms; */
/*   sum = 0; */
/*   maxFeatures = (int*)malloc(sizeof(int)*m->mode); */

/*   for(i = 0; i < testSet->n; i++){ */

/*     ////////////////////////////////////////////////// Computing maxFeatures ////////////////////////////////////////////////////// */
/*     for(j = 0; j < m->mode; j++){ */
/*       maxFeatures[j] = (testSet->features)[i]; */
/*       for(k = 0; k < (testSet->features)[i]; k++){ */
/* 	for(l = 0; l < (m->mWidth)[j]; l++){ */
/* 	  if(k + l > (testSet->features)[i] - (m->mWidth)[j]){ */
/* 	    maxFeatures[j]--; */
/* 	    break; */
/* 	  } */
/* 	  if((testSet->data)[i][k + l] > 3){ */
/* 	    maxFeatures[j]--; */
/* 	    break; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */

/*     max = 0; */
/*     for(j = 0; j < (testSet->features)[i]; j++){ */

/*       //////////////////////////////////////// mode count ///////////////////////////////////// */
/*       modeCount = m->mode; */
/*       for(k = 0; k < m->mode; k++){ */
/* 	for(l = 0; l < (m->mWidth)[k]; l++){ */
/* 	  if(j + l > (testSet->features)[i] - (m->mWidth)[k]){ */
/* 	    modeCount--; */
/* 	    break; */
/* 	  } */
/* 	  if((testSet->data)[i][j + l] > 3){ */
/* 	    modeCount--; */
/* 	    break; */
/* 	  } */
/* 	} */
/*       } */
/*       ///////////////////////////////////////////////////////////////////////////////////////// */
/*       for(k = 0; k < m->mode; k++){ */
/* 	//	s2 = ((m->t)[k] + m->featureValues*m->alpha) / (m->n + m->mode*m->featureValues*m->alpha); */
/* 	//	s2 = 1; */
/* 	s2 = 0; */
/* 	ms = (m->motifs)[k].motif; */
/* 	for(l = 0; l < (m->mWidth)[k]; l++){ */
/* 	  if(j + l > (testSet->features)[i] - (m->mWidth)[k]){ */
/* 	    s2 = 0; */
/* 	    break; */
/* 	  } */
/* 	  if((testSet->data)[i][j + l] > 3){ */
/* 	    s2 = 0; */
/* 	    break; */
/* 	  } */
/* 	  //	  s2 = s2 * ((ms->modeMotifCount)[(testSet->data)[i][j + l]] + m->alpha) / ((m->t)[k] + m->featureValues*m->alpha); */
/* 	  s2 = s2 + log((ms->modeMotifCount)[(testSet->data)[i][j + l]] + m->alpha) - log((m->t)[k] + m->featureValues*m->alpha); */
/* 	  //	  if(background[i][j + l] > 0.00001) s2 = s2 / background[i][j + l]; */
/* 	  if(background[i][j + l] > 0.00001) s2 = s2 - log(background[i][j + l]); */
/* 	  ms = ms->next; */
/* 	} */

/* 	if(s2 == 0) continue; */

/* 	if(modeCount == 0 && s2 != 0) printf("ERROR!!!\n"); */
/* 	if(maxFeatures[k] == 0 && s2 != 0) printf("ERROR1!!!\n"); */

/* 	if(modeCount != 0) s2 = s2 - log((double)modeCount); */
/* 	if(maxFeatures[k] != 0) s2 = s2 - log((double)maxFeatures[k]); */

/* 	if(s2 > max) max = s2; */
/*       } */
/*     } */
/*     sum = sum + max; */
/*     //    sum = sum + log(max); */
/*   } */
/*   free(maxFeatures); */
/*   return sum; */
/* } */
