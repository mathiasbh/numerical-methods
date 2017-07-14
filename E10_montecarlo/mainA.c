#include <math.h>
#include <stdio.h>
#include "header.h"
#define M_PI 3.14159265358979323846264338


double ftest1(double* x) {return x[0]+x[1];}
double ftest2(double* x) {return x[0]*x[1]*x[2]*x[3];}

double f(double *x) {
  double r = 1.0-cos(x[0])*cos(x[1])*cos(x[2]);
  return 1.0/r/M_PI/M_PI/M_PI;
}

int main() {
  int size = 2, N = 1e4;
  double error, result, expected;

  fprintf(stderr, "#################### Exercise A ####################\n");
  double a1[] = {0.0, 0.0}, b1[] = {2.0, 2.0};
  expected = 8;
  plainmc(size,ftest1,a1,b1,N,&result,&error);
  fprintf(stderr," ---- Integral x+y, x,y=[%lg,%lg]\n", a1[0], b1[0]);
  fprintf(stderr,"Result = %lg, error = %lg, points N = %i, expected = %lg\n\n", result ,error, N, expected);

  size = 4; N = 1e6; expected = 20.25;
  double a2[] = {0.0, 0.0, 0.0, 0.0}, b2[] = {2.0, 2.0, 1.5, 3.0};
  plainmc(size,ftest2,a2,b2,N,&result,&error);
  fprintf(stderr," ---- Integral x*y*z*w, x,y=[%lg,%lg], z=[%lg,%lg], w=[%lg,%lg]\n", a2[0],b2[0],a2[2],b2[2],a2[3],b2[3]);
  fprintf(stderr,"Result = %lg, error = %lg, points N = %i, expected = %lg\n\n", result ,error, N, expected);

  /* Test if error goes as 1/sqrt(N) */
  size = 2; N = 1e4;
  for(int i=2; i<200; i++) {
    plainmc(size,ftest1,a1,b1,i,&result,&error);
    printf("%i %lg\n", i, error);
  }

  /* Calculate given function, x,y,z=[0,PI] - should be fres: */
  double fres = 1.3932039296856768591842462603255;
  double a3[] = {0.0, 0.0, 0.0}, b3[] = {M_PI, M_PI, M_PI};
  size = 3;
  N = 1e5; // use 1e8 to get first three digits correctly
  plainmc(size,f,a3,b3,N,&result,&error);
  fprintf(stderr," ---- Integral pi^(-3)/(1-cos(x)cos(y)cos(z)), x,y,z=[0,PI]\n");
  fprintf(stderr,"Result = %0.7lg, error = %lg, expected = %0.7lg, points N = %i\n\n", result ,error, fres, N);

  return 0;
}
