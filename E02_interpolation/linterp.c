/* 
  Linear interpolation through x[i], y[i] at x = z
  int n: number of data 
  double x: array
  double y: array
  double z: number

  Name: Mathias Bojsen-Hansen
  DATE: 12-04-2016
*/

#include <assert.h>

double linterp(int n, double *x, double *y, double z) {
  /* Require more than one data point, z between first and last datapoint */
  assert(n>1);
  assert(z>=x[0] && z<=x[n-1]);

  int i = binary_search(n, x, z);

  return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i]) * (z - x[i]);
}
