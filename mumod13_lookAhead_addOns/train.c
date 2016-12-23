#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "messages.h"
#include "linkedListOps.h"
#include "modelops.h"
#include "motifops.h"
#include "traindata.h"
#include "crossValidation.h"
int main(){

  dataSet *ds, *dsT, *dsL;
  model *m;
  int *labels, *startPos, mode, i, j, k, k1, *lpc, *spc, minWidth, maxWidth;
  unsigned int seed;
  float alpha;
  double l1, zoops;
  double **background;
  int **featureCounts;
  FILE *fp, *fp3;
  char *fstr;
  trainOut *l;
  int kfold;
  int *posL, *mWidth, bestMode;
  trainOut *savedTo;
  double cvL, bestCVL, tmpCVL, tmpL, like;
  char *strName;
  mode = 8;
  fstr = "likes.txt";
  float lambda;
  /* int mWidth[] = {8, 8, 8}; */
  alpha = 1;
  zoops = 0;
  lambda = 0;
  bestCVL = 0;
  minWidth = 6;
  maxWidth = 15;
  kfold = 5;
  //  ds = getData("new_data_1.txt", "test.txt", 1, 1);
  ds = getData("../mumodData/GSM602326_dmel_summit_201bp.fa", "test.txt", 1, 1);
  /* for(i = 0; i < ds->n; i++){ */
  /*   for(j = 0; j < (ds->features)[i]; j++) */
  /*     printf("%d\t", (ds->lookAhead)[i][j]); */
  /*   printf("\n"); */
  /* } */
  /* exit(0); */
  //  ds = getData("../../mumodData/GSM602326_dmel_summit_201bp.fa", "test.txt", 1);
  //  ds = getData("../../mumodData/GSM602330_dsim_summit_201bp.fa", "test.txt", 1);
  //  ds = getData("../../mumodData/GSM602333_dyak_summit_201bp.fa", "test.txt", 1);
  //  ds = getData("../../mumodData/GSM602335_dpse_summit_201bp.fa", "test.txt", 1);

  //  featureCounts = getFeatureCounts(ds, mode, 8);
  background = getBackground1(ds, 2);

  /* printf("Counts:\n"); */
  /* for(i = 0; i < ds->n; i++){ */
  /*   for(j = 0; j < mode; j++){ */
  /*     printf("%d\t", featureCounts[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* exit(0); */
  l1 = 0;

  mWidth = (int*)malloc(sizeof(int)*mode);
  for(i = 0; i < mode; i++)
    mWidth[i] = 8;

  for(k = 0; k < 1; k++){
    printf("Trial: %d\n", k);
    seed = k + 1;
    l = trainData(ds, mode, alpha, lambda, zoops, seed, background, mWidth, minWidth, fstr);
    printf("Likelihood: %lf\n", l->likelihood);
    printf("Widths:\n");
    for(i = 0; i < mode; i++) printf("%d\t", (l->motifWidth)[i]);
    printf("\n");
    if(l1 == 0 || l->likelihood > l1){
      l1 = l->likelihood;
      savedTo = l;
    }
  }
  fp = fopen("learnedLabels.txt", "w");
  for(i = 0; i < ds->n; i++) fprintf(fp, "%d\n", (savedTo->labels)[i]);
  fclose(fp);
  fp = fopen("learnedStart.txt", "w");
  for(i = 0; i < ds->n; i++) fprintf(fp, "%d\n", (savedTo->startPos)[i]);
  fclose(fp);
  fp = fopen("learnedWidth.txt", "w");
  for(i = 0; i < mode; i++) fprintf(fp, "%d\n", (savedTo->motifWidth)[i]);
  fclose(fp);

  exit(0);
  posL = posList(ds->n);

  for(i = 2; i < 11; i++){
    printf("MODEL: %d modes\n", i);
    mWidth = (int*)malloc(sizeof(int)*i);
    cvL = 0;
    like = 0;
    for(j = 0; j < i; j++) mWidth[j] = 8;
    for(j = 0; j < kfold; j++){
      printf("Fold: %d\n", j);
      dsT = getTrainSubset(ds, j, kfold, posL);
      /* featureCounts = getFeatureCounts(dsT, minWidth, maxWidth); */
      /* for(k = 0; k < dsT->n; k++){ */
      /* 	for(k1 = 0; k1 < maxWidth - minWidth + 1; k1++){ */
      /* 	  printf("%d\t", featureCounts[k][k1]); */
      /* 	} */
      /* 	printf("\n"); */
      /* } */
      /* exit(0); */

      background = getBackground1(dsT, 2);
      l1 = 0;
      for(k = 0; k < 5; k++){
	printf("Trial: %d\n", k);
	seed = k + 1;
	l = trainData(dsT, i, alpha, lambda, zoops, seed, background, mWidth, minWidth, fstr);
	//	printf("HERE!!\n");
	if(k == 0){
	  savedTo = l;
	}
	else if(l->likelihood > l1){
	  l1 = l->likelihood;
	  freeTo(savedTo);
	  savedTo = l;
	}
	else{
	  freeTo(l);
	}
      }

      for(k = 0; k < dsT->n; k++) free(background[k]);
      free(background);
      /* for(k = 0; k < dsT->n; k++) free(featureCounts[k]); */
      /* free(featureCounts); */
      dsL = getTestSubset(ds, j, kfold, posL);
      //      dsL = getTrainSubset(ds, j, kfold, posL);
      background = getBackground1(dsL, 2);
      ///////////////////////////////////////////////// STR COPY ////////////////////////////////////////////////////

      strName = malloc(sizeof(char)*(9 + 6 + 4 + (floor(log10(abs(i)))) + 1 + 1));
      strName[0] = '\0';
      sprintf(strName, "mode%dfold%dlabels.txt", i, j);
      printf("%s\n", strName);
      fp = fopen(strName, "w");
      for(k = 0; k < dsT->n; k++) fprintf(fp, "%d\n", (savedTo->labels)[k]);
      fclose(fp);
      free(strName);

      strName = malloc(sizeof(char)*(9 + 5 + 4 + (floor(log10(abs(i)))) + 1 + 1));
      strName[0] = '\0';
      sprintf(strName, "mode%dfold%dstart.txt", i, j);
      printf("%s\n", strName);
      fp = fopen(strName, "w");
      for(k = 0; k < dsT->n; k++) fprintf(fp, "%d\n", (savedTo->startPos)[k]);
      fclose(fp);
      free(strName);

      strName = malloc(sizeof(char)*(9 + 5 + 4 + (floor(log10(abs(i)))) + 1 + 1));
      strName[0] = '\0';
      sprintf(strName, "mode%dfold%dwidth.txt", i, j);
      printf("%s\n", strName);
      fp = fopen(strName, "w");
      for(k = 0; k < i; k++) fprintf(fp, "%d\n", (savedTo->motifWidth)[k]);
      fclose(fp);
      free(strName);
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      /* m = createModel(i, dsT, savedTo->labels, savedTo->startPos, alpha, lambda, zoops, savedTo->motifWidth); */

      /* printf("Model mode counts:\n"); */
      /* for(k = 0; k < m->mode; k++) printf("%d\t", (m->t)[k]); */
      /* printf("\n"); */
      /* printf("Learnt model motif widths:\n"); */
      /* for(k = 0; k < m->mode; k++) printf("%d\t", (m->mWidth)[k]); */
      /* printf("\n"); */

      /* tmpCVL = cvLikelihood(dsL, m, background); */
      /* tmpL = savedTo->likelihood; */
      /* printf("Likelihood: %lf\n", savedTo->likelihood); */
      /* printf("CVL out: %lf\n", tmpCVL); */
      /* cvL = cvL + tmpCVL; */
      /* like = like + tmpL; */
      /* //      printf("CVL: %lf\n", cvL); */
      /* freeModel(m); */
      freeData(dsL);
      freeData(dsT);
      freeTo(savedTo);
    }
    cvL = cvL/5;
    like = like/5;
    printf("Average: %lf\n", cvL);
    printf("Average Likelihood: %lf\n", like);
    if(bestCVL == 0 || cvL > bestCVL){
      bestCVL = cvL;
      bestMode = i;
    }
    free(mWidth);
  }
  printf("Best MODEL: %d modes\n", bestMode);
  freeData(ds);


  /* background = getBackground1(ds, 2); */
  /* l1 = 0; */
  /* for(i = 0; i < 5; i++){ */
  /*   seed = i + 1; */
  /*   printf("Trial: %d\n", i + 1); */
  /*   printf("Zoops: %lf\n", zoops); */
  /*   if(i == 0) fstr = "likes1.txt"; */
  /*   if(i == 1) fstr = "likes2.txt"; */
  /*   if(i == 2) fstr = "likes3.txt"; */
  /*   if(i == 3) fstr = "likes4.txt"; */
  /*   if(i == 4) fstr = "likes5.txt"; */
  /*   l = trainData(ds, mode, alpha, lambda, zoops, seed, background, mWidth, minWidth, fstr); */
  /*   if(i == 0 || l->likelihood > l1){ */
  /*     printf("Saving trial: %d\n", i + 1); */
  /*     fp = fopen("learnedLabels.txt", "w"); */
  /*     for(j = 0; j < ds->n; j++){ */
  /* 	fprintf(fp, "%d\n", (l->labels)[j] + 1); */
  /*     } */
  /*     fclose(fp); */

  /*     fp = fopen("learnedStart.txt", "w"); */
  /*     for(j = 0; j < ds->n; j++){ */
  /*   	fprintf(fp, "%d\n", (l->startPos)[j]); */
  /*     } */
  /*     fclose(fp); */

  /*     fp = fopen("learnedWidth.txt", "w"); */
  /*     for(j = 0; j < mode; j++) */
  /*   	fprintf(fp, "%d\n", (l->motifWidth)[j]); */
  /*     fclose(fp); */

  /*     l1 = l->likelihood; */
  /*   } */

  /*   freeTo(l); */
  /* } */

  /* for(i = 0; i < ds->n; i++) free(background[i]); */
  /* free(background); */
  /* freeData(ds); */
  return 0;
}
