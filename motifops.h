#ifndef _motifops_h
#define _motifops_h

double motifAddLeftScore(dataSet*, model*, double**, int*, int*, int);
double motifRemoveLeftScore(dataSet*, model*, double**, int*, int*, int);
double motifAddRightScore(dataSet*, model*, double**, int*, int*, int);
double motifRemoveRightScore(dataSet*, model*, double**, int*, int*, int);
void motifLeftIncrease(dataSet*, model*, double**, int*, int*, int);
void motifLeftDecrease(dataSet*, model*, double**, int*, int*, int);
void motifRightIncrease(dataSet*, model*, double**, int*, int*, int);
void motifRightDecrease(dataSet*, model*, double**, int*, int*, int);

#endif
