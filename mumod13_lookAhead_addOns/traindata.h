#ifndef _traindata_h
#define _traindata_h

typedef struct traindatastruct{
  model *m;
  dataSet *ds;
  int *labels;
  int *startPos;
  int times;
  double likelihood;
  unsigned int seed;
  char *filename;
}trainDataStruct;

typedef struct queuestruct{
  double y;
  double xy;
}queueStruct;

double linearSlope(queueStruct*, double, double, double, double, int);
double corrCoef(int*, double*, int);
double avgCorrCoef(model*, dataSet*, int*, int*, double**, int**);
void enqueue(queueStruct*, int*, int*, int*, int, double, double*, double*);
void dequeue(queueStruct*, int*, int*, int, double*, double*);
void addRemoveDataPoint(model*, dataSet*, int*, int*, int, int);
int updateBestModel(dataSet*, model*, double**, int*, int*, int);
int bestPosteriorParameters(dataSet*, model*, double*, int*, int*, int, int*, int, int);
double calculateFullPosterior(dataSet*, model*, double**, int*, int*, int**);
double calculateLikelihood(model*, dataSet*, int*, int*, double**, int**, int);
double calculateLikelihoodMode(model*, dataSet*, int*, int*, double**, int**, int, int);
double calculateLikelihoodFaster(model*);
double calculateLikelihoodSingleMode(model*, int);
int sampleNewLabel(model*, dataSet*, int*, int, unsigned int*, int, double*, int*, int);
int sampleStartPosFast(model*, dataSet*, int, int, int, unsigned int*, int, double*);
int sampleStartPos(model*, dataSet*, int, int, unsigned int*, int, double*);
int sampleMotifWidthLeft(dataSet*, model*, double**, int*, int*, int, int, unsigned int*);
int sampleMotifWidthRight(dataSet*, model*, double**, int*, int*, int, int, unsigned int*);
double EMLike(model*, dataSet*, int*, int*, double**);
trainOut* trainData(dataSet*, int, float, float, double, unsigned int, double**, int*, int, char*);
double likelihoodXi(model*, dataSet*, int, int, int, double*);
double likelihoodXiFaster(model*, dataSet*, int, int, int, double*);
int* getStartIndices(model*, dataSet*, int, int, unsigned int*, int, double*);
#endif