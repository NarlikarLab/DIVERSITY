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

#ifndef _modelops_h
#define _modelops_h

typedef struct dataset{
  int **data;
  int *features; 
  int **lookAhead;
  int featureValues;
  int n;
  int tu;
}dataSet;

typedef struct modelstruct{
  int mode;
  int featureValues;
  int n;
  motifContainer *motifs;
  int *t;
  int *mWidth;
  float alpha;
  float lambda;
  double zoops;
}model;

typedef struct trainout{
  int *labels;
  int *startPos;
  int *motifWidth;
  double likelihood;
}trainOut;

dataSet* getData(char*, char*, int, int);
double** getBackground(dataSet*, int);
int maxWithoutN(int*, int, int);
int getIndex(int*, int, int, int);
void freeData(dataSet*);
model* createModel(int, dataSet*, int*, int*, float, float, double, int*);
model* initializeModel(int, int, int*);
void copyModel(model*, model*);
void freeModel(model*);
void freeTo(trainOut*);
void getMotifCount(model*, dataSet*, int*, int*);
int* sampleArrayVals(int, int, int, unsigned int*);
int* sampleStart(int, int, int, int, unsigned int*);
void logToExp(double*, int);
int sample(double*, int, double);
void bubbleSort(double*, int*, int);
int arrayMax(double*, int);
void copyLabels(int*, int*, int);
int motifWithN(int*, int, int);
void initializeLabelStartPos(dataSet*, int*, int*, int, int*, unsigned int*);
int* nRandomPos(int, int, int, unsigned int*);
double* normalize(int*, int);
int min(int, int);
int** getFeatureCounts(dataSet*, int, int);
void modeFeatureCount(dataSet*, int**, int, int);
int** lookForN(dataSet*);
#endif
