#include <math.h>
#include <stdio.h>
#include "header.h"
#define M_PI 3.14159265358979323846264338

int calls = 0;

double ftest(double* x) {calls++; return x[0]*x[1]*x[1]*x[1];}

double ftest2(double *x) {
  calls++;
  return x[0]*x[1]/(1.0+x[0]*x[1]*x[1]);
}

int main() {
  int size = 2;
  double error, result, acc=1e-5, eps=acc, expected;
  double a[] = {0.0, 0.0}, b[] = {2.5, 2.5};
  double tol = 5e-3;
  
  fprintf(stderr, "#################### Exercise B ####################\n");
  double Q = intOpen2d(size,ftest,a,b,acc,eps,&error); expected = 30.5176; error = fabs(expected - Q);
  fprintf(stderr," ---- Integral (recursive 1D adaptive) x*y^3, x,y=[%lg,%lg]\n", a[0],b[0]);
  fprintf(stderr,"Result = %lg, error = %lg, expected = %lg, calls=%i\n\n",Q,error,expected,calls);

  calls = 0;
  mctol(size,ftest,a,b,tol,&result,&error); 
  fprintf(stderr," ---- Integral (MC using tolerance) x*y^3, x,y=[%lg,%lg]\n", a[0],b[0]);
  fprintf(stderr,"Result = %lg, error = %lg, tolerance = %lg, expected = %lg, calls=%i\n\n",result,error,tol,expected,calls);

  double a1[] = {0.0, 0.0}, b1[] = {2.0, 2.0};
  error =3.1415; calls=0;
  Q = intOpen2d(size,ftest2,a1,b1,acc,eps,&error); expected = 9/8.0 * log(9) - 1; error = fabs(expected - Q);
  fprintf(stderr," ---- Integral (recursive 1D adaptive) xy/(1+x*y^2), x,y=[%lg,%lg]\n", a1[0],b1[0]);
  fprintf(stderr,"Result = %lg, error = %lg, expected = %lg, calls=%i\n\n",Q,error,expected,calls);

  calls=0;
  mctol(size,ftest2,a1,b1,tol,&result,&error);
  fprintf(stderr," ---- Integral (MC using tolerance) xy/(1+x*y^2), x,y=[%lg,%lg]\n", a1[0],b1[0]);
  fprintf(stderr,"Result = %lg, error = %lg, tolerance = %lg, expected = %lg, calls=%i\n\n",result,error,tol,expected,calls);
  return 0;
}
