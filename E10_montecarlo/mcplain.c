#include <stdio.h>
#include <math.h>

void samplingRand(int size, double* a, double* b, double *x);
  
void plainmc(int size, double f(double* x), double* a, double* b, int N, double* result, double* error) {
  double x[size], fx, sumfx = 0, sumProd = 0, volume = 1;
  for(int i=0; i<size; i++) {
    volume *= (b[i] - a[i]);
  }
  for(int i=0; i<N; i++) {
    samplingRand(size, a, b, x);
    fx = f(x);
    sumfx += fx;
    sumProd += fx*fx;
  }
  double meanVal = sumfx/N;
  double var = sumProd/N - meanVal*meanVal; // eqn 4
  *result = meanVal*volume;
  *error = sqrt(var/(double)N)*volume; // eqn 3
}
