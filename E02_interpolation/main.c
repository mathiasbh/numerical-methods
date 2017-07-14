/*
  Perform interpolation using
  Linear method
  Quadratic spline (+derivative and integral)
  Cubic spline (+derivative and integral)
*/

#include <stdio.h>
#include <math.h>
#include "cspline.h"
#include "qspline.h"
#define M_PI 3.14159265358979323846 

double linterp(int n, double *x, double *y, double z);


int main() {
  /* Allocate memory to arrays */
  int n = 10;
  double x[n], y[n];

  /* Function to interpolate */
  for(int i=0; i<n; i++) {
    x[i] = -4.6+i;
    y[i] = exp(-x[i]*x[i]);
    printf("%g\t%g\n", x[i], y[i]);
  } printf("\n\n");

  double step = (x[n-1]-x[0])/100.0;
  /* Exact solutions, function, derivative, integral */
  for(double z=x[0]; z<=x[n-1]; z+=step) {
    printf("%lg \t %lg \t %lg \t %lg\n", z, exp(-z*z), -2*z*exp(-z*z), sqrt(M_PI)*(erf(z)-erf(x[0]))/2);
  } printf("\n\n");

  /* Interpolate data using linear interpolation */
  for(double z=x[0]; z<=x[n-1]; z+=step) {
    printf("%lg\t%lg\n", z, linterp(n,x,y,z));
  } printf("\n\n");

  /* Interpolate data using quadratic interpolation */
  qspline* qs = qspline_alloc(n, x, y);
  for(double z=x[0]; z<=x[n-1]; z+=step) {
    printf("%lg\t%lg\t%lg\t%lg\n", z, qspline_evaluate(qs,z), qspline_derivative(qs,z), qspline_integral(qs,z));
  } printf("\n\n");

  qspline_free(qs);

  /* Interpolate data using cubic interpolation */
  cspline* cs = cspline_alloc(n, x, y);
  for(double z=x[0]; z<=x[n-1]; z+=step) {
    printf("%lg\t%lg\t%lg\t%lg\n", z, cspline_evaluate(cs,z), cspline_derivative(cs,z), cspline_integral(cs,z));
  } printf("\n\n");

  cspline_free(cs);
  
  return 0;
}
