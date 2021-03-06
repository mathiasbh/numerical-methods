#include <stdio.h>
#include <math.h>
#include "header.h"
#include "header1.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>


int main() {
  int calls;
  double a, b, acc, eps, Q, err, expected;
  double fpwr2  (double x) {calls++; return x*x;};
  double fmbd   (double x) {calls++; return x*x*exp(-x*x);}
  double finv   (double x) {calls++; return 1.0/sqrt(1.0+x);}
  double f      (double x) {calls++; return 4.0*sqrt(1.0-(1.0-x)*(1.0-x));}
  double fpi0inf(double x) {calls++; return 1.0/(x*x+2*2);} // 0 to inf = pi/4=0.785398
  double fexpx2 (double x) {calls++; return exp(-2*x*x);} // -inf to inf = sqrt(pi/2) = 1.25331
  double finvx22(double x) {calls++; return 1.0/(2.0+x*x);} // -inf to 0 = pi/(2sqrt2) = 1.1107
  void fode(double x, gsl_vector*y, gsl_vector*dydt) {gsl_vector_set(dydt,0,4.0*sqrt(1.0-(1.0-x)*(1.0-x)));}
  void fode1(double x, gsl_vector*y, gsl_vector*dydt) {gsl_vector_set(dydt,0,x*x);}
  
  printf("##################Exercise A##################\n");
  a=-1; b=5; acc=1e-4; eps=1e-4; calls=0; expected = 42;
  Q=intOpen(fpwr2,a,b,acc,eps,&err);
  printf("--- Integral f=x^2 (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i, expected = %lg\n\n\n", Q, err, calls, expected);
  
  calls=0;
  Q=intClosed(fpwr2,a,b,acc,eps,&err);
  printf("--- Integral f=x^2 (a,b)=(%lg,%lg): closed quadrature ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i, expected = %lg\n\n\n", Q, err, calls, expected);

  
  a=-10; b=8; calls=0; expected = 0.886227;
  Q=intOpen(fmbd,a,b,acc,eps,&err);
  printf("--- Integral f=x^2*e^{-x^2} (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i, expected = %lg\n\n\n", Q, err, calls, expected);

  calls=0;
  Q=intClosed(fmbd,a,b,acc,eps,&err);
  printf("--- Integral f=x^2*e^{-x^2} (a,b)=(%lg,%lg): closed quadrature ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i, expected = %lg\n\n\n", Q, err, calls, expected);

  
  a=0; b=1; acc=1e-15; eps=acc; calls=0;
  Q=intOpen(f,a,b,acc,eps,&err);
  printf("--- Integral f=4*sq(1-(1-x)^2) (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%.17lg, err=%lg, calls=%i, acc=eps=%lg\n\n\n", Q, err, calls,acc);

  calls=0;
  Q=intClosed(f,a,b,acc,eps,&err);
  printf("--- Integral f=4*sq(1-(1-x)^2) (a,b)=(%lg,%lg): closed quadrature ---\n",a,b);
  printf("Q=%.17lg, err=%lg, calls=%i, acc=eps=%lg\n\n\n", Q, err, calls,acc);

  double fgsl (double x, void *params) {calls++;return 4.0*sqrt(1.0-(1.0-x)*(1.0-x));}
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  gsl_function F = {.function = fgsl, .params=nan};
  calls=0; acc = 1e-13; eps=1e-16;
  gsl_integration_qags (&F, a, b, acc, eps, 10000, w, &result, &error); 
  printf("--- Integral f=4*sq(1-(1-x)^2): GSL QAGS ---\n");
  printf("Q=%.17lg, err=%lg, calls=%i\n\n\n", result, error,calls);
  gsl_integration_workspace_free (w);
  printf("gsl is much faster but will not go to as high precision due to roundoff errors,\n");
  printf(" however, we only get one more digit of pi with 6000 times more calls...\n\n\n");


  printf("##################Exercise B##################\n");
  /* B: Infinity integrals */
  acc=1e-4; eps=1e-4;

  calls=0; a=0; b=INFINITY;
  Q=intOpen(fpi0inf,a,b,acc,eps,&err);
  printf("--- Integral f=1.0/(x^2+2^2) (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%.7lg (should be 0.785398), err=%lg, calls=%i\n\n\n", Q, err, calls);


  calls=0; a=-INFINITY; b=INFINITY;
  Q=intOpen(fexpx2,a,b,acc,eps,&err);
  printf("--- Integral f=exp(-2*x^2) (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%.7lg (should be 1.25331), err=%lg, calls=%i\n\n\n", Q, err, calls);


  calls=0; a=-INFINITY; b=0;
  Q=intOpen(finvx22,a,b,acc,eps,&err);
  printf("--- Integral f=1.0/(2.0+x*x) (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%.7lg (should be 1.110721), err=%lg, calls=%i\n\n\n", Q, err, calls);

  
  printf("##################Exercise C##################\n");
  /* C: Infinity integrals */
  a=-1; b=1; acc=1e-6; eps=1e-6; calls=0;
  Q = clensshawcurtis(finv,a,b,acc,eps,&err);
  printf("--- Integral f=x^2 (a,b)=(%lg,%lg): clensshawcurtis ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i\n\n\n", Q, err, calls);


  a=-1; b=1; calls=0;
  Q = intOpen(finv,a,b,acc,eps,&err);
  printf("--- Integral f=x^2 (a,b)=(%lg,%lg): without clensshawcurtis ---\n",a,b);
  printf("Q=%lg, err=%lg, calls=%i\n\n\n", Q, err, calls);

  /* C: using ODE, compare with open and closed adaptive */
  int n = 1; double t=0; a=0; b=1; double h=0.1; acc=1e-8; eps=acc;
  gsl_vector* y = gsl_vector_alloc(n); gsl_vector_set(y, 0, a); 
  calls = driver(&t,b,&h,y,acc,eps,rkstep5,fode);
  printf("--- Integral f=4*sq(1-(1-x)^2) (a,b)=(%lg,%lg): using ODE routines ---\n",a,b);
  printf("Q=%.7lg, calls=%i, acc=eps=%lg\n\n\n", gsl_vector_get(y,0), calls,acc);
  gsl_vector_free(y);

  calls=0;
  Q=intOpen(f,a,b,acc,eps,&err);
  printf("--- Integral f=4*sq(1-(1-x)^2) (a,b)=(%lg,%lg): open quadrature ---\n",a,b);
  printf("Q=%.7lg, err=%lg, calls=%i, acc=eps=%lg\n\n\n", Q, err, calls,acc);

  calls=0;
  Q=intClosed(f,a,b,acc,eps,&err);
  printf("--- Integral f=4*sq(1-(1-x)^2) (a,b)=(%lg,%lg): closed quadrature ---\n",a,b);
  printf("Q=%.7lg, err=%lg, calls=%i, acc=eps=%lg\n\n\n", Q, err, calls,acc);
  return 0;
}
