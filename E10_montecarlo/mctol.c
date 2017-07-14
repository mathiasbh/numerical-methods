#include <stdio.h>
#include <math.h>

void samplingRand(int size, double* a, double* b, double *x);

void mctol(int size, double f(double* x), double* a, double* b, double tol, double* result, double* error) {
  int N = 1000, Ntotal = 0;
  double x[size], fx, sumfx = 0, sumProd = 0, volume = 1, meanVal, var, sigma;
  for(int i=0; i<size; i++) {
    volume *= (b[i] - a[i]);
  }
  while(1) {
    Ntotal+=N;
    for(int i=0; i<N; i++) {
      samplingRand(size, a, b, x);
      fx = f(x);
      sumfx += fx;
      sumProd += fx*fx;
    }
    meanVal = sumfx/Ntotal;
    var = sumProd/Ntotal - meanVal*meanVal; // eqn 4
    sigma = sqrt(var/(double)Ntotal);
    if(sigma*volume < tol) break;
  }
  *result = meanVal*volume;
  *error = sigma*volume; // eqn 3
}
