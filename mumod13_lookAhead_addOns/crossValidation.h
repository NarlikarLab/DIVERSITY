#ifndef _crossValidation_h
#define _crossValidation_h

void nRandom1(int*, int, int, unsigned int*);
int* posList(int);
dataSet* getTrainSubset(dataSet*, int, int, int*);
dataSet* getTestSubset(dataSet*, int, int, int*);
double chooseBestStart(dataSet*, model*, double**, int, int);
double cvLikelihood(dataSet*, model*, double**);
#endif
