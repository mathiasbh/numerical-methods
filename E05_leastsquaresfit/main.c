
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "header.h"
#define RND ((double)rand()/RAND_MAX)

double funcA(int i, double x) {
  switch(i) {
  case 0: return 1.0/x; break;
  case 1: return 1.0;   break;
  case 2: return x;     break;
  default: {fprintf(stderr, "Function A: wrong input i: %d\n", i); return NAN;}
  }
}

double funcB(int i, double x) {
  switch(i) {
  case 0: return 1.0;       break;
  case 1: return 1.0*x;     break;
  case 2: return x*sin(x);  break;
  default: {fprintf(stderr, "Function B: wrong input i: %d\n", i); return NAN;}
  }
}

double experimental_funcB(double x) {return 3 + 0.4*x + 0.5*sin(x)*x;}

int main() {

  /* Exercise A */
  int n = 10; /* Number of datapoints */
  double xdata[]  = {0.100, 0.145, 0.211, 0.307, 0.447, 0.649, 0.944, 1.372, 1.995, 2.900};
  double ydata[]  = {12.644, 9.235, 7.377, 6.460, 5.555, 5.896, 5.673, 6.964, 8.896, 11.355};
  double dydata[] = {0.858, 0.359, 0.505, 0.403, 0.683, 0.605, 0.856, 0.351, 1.083, 1.002};

  gsl_vector* x  = gsl_vector_alloc(n);
  gsl_vector* y  = gsl_vector_alloc(n);
  gsl_vector* dy = gsl_vector_alloc(n);

  for(int i=0; i<n; i++) {
    gsl_vector_set(x,  i, xdata[i]);
    gsl_vector_set(y,  i, ydata[i]);
    gsl_vector_set(dy, i, dydata[i]);
    printf("%lg\t%lg\t%lg\n", xdata[i], ydata[i], dydata[i]);
  } printf("\n\n");

  int m = 3; /* Three parameters */
  gsl_vector* c  = gsl_vector_alloc(m);
  gsl_vector* dc = gsl_vector_alloc(m);
  gsl_matrix* S  = gsl_matrix_alloc(m,m);

  lsfit(x,y,dy,m,funcA,c,S); /* Returns c in Ac = b and S (uncertainty) */

  /* Uncertainty, Delta c */
  fprintf(stderr, "Exercise A fitting coefficients\n");
  for(int i=0; i<m; i++) gsl_vector_set(dc, i, sqrt(gsl_matrix_get(S,i,i)));

  /* Print result */
  fprintf(stderr,"f(x) = a1/x + a2 + a3*x\n[a1,a2,a3] = \n");
  for(int i=0; i<m; i++) fprintf(stderr,"%7.5lg +/- %7.5lg\n", gsl_vector_get(c,i), gsl_vector_get(dc,i)); 
  fprintf(stderr,"\n");

  /* Print function for plotting */
  double a1 = gsl_vector_get(c,0), a2 = gsl_vector_get(c,1), a3 = gsl_vector_get(c,2);
  double x0=0.07, x1=3.0, dx=0.05;
  for(double x=x0; x<x1; x+=dx) {

    double fx = a1/x  + a2  + a3*x;
    double fxU= (a1+gsl_vector_get(dc,0))/x + (a2+gsl_vector_get(dc,1)) + (a3+gsl_vector_get(dc,2))*x;
    double fxL= (a1-gsl_vector_get(dc,0))/x + (a2-gsl_vector_get(dc,1)) + (a3-gsl_vector_get(dc,2))*x;
    printf("%lg\t%lg\t%lg\t%lg\n", x, fx, fxU, fxL);
  } printf("\n\n");




  /* Exercise B */

  /* n = 10, already defined */
  double coeff1=0, coeff2=16;
  for(int i=0; i<n; i++) {
    double xi = coeff1 + (coeff2-coeff1)*i/(n-1);
    gsl_vector_set(x,  i, xi);
    gsl_vector_set(y,  i, experimental_funcB(xi) + 0.4*RND);
    gsl_vector_set(dy, i, 1.3*RND);
    printf("%lg\t%lg\t%lg\n", gsl_vector_get(x,i), gsl_vector_get(y,i), gsl_vector_get(dy,i));
  } printf("\n\n");

  /* m = 3, already defined */
  /* gsl_vector* c  = gsl_vector_alloc(m); */
  /* gsl_vector* dc = gsl_vector_alloc(m); */
  /* gsl_matrix* S  = gsl_matrix_alloc(m,m); */

  lsfit(x,y,dy,m,funcB,c,S); /* Returns c in Ac = b and S (uncertainty) */

  /* Uncertainty, Delta c */
  fprintf(stderr, "Exercise B fitting coefficients\n");
  for(int i=0; i<m; i++) gsl_vector_set(dc, i, sqrt(gsl_matrix_get(S,i,i)));

  /* Print result */
  fprintf(stderr,"f(x) = a1 + a2*x + a3*x*sin(x)\n[a1,a2,a3] = \n");
  for(int i=0; i<m; i++) fprintf(stderr,"%7.5lg +/- %7.5lg\n", gsl_vector_get(c,i), gsl_vector_get(dc,i)); 
  fprintf(stderr, "Input: a1=3.0, a2=0.4, a3=0.5\n"); fprintf(stderr,"\n");

  /* Construct numbers for plotting */
  a1 = gsl_vector_get(c,0), a2 = gsl_vector_get(c,1), a3 = gsl_vector_get(c,2);
  dx=0.05;
  for(double x=coeff1; x<coeff2; x+=dx) {
    double fx = a1  + a2*x  + a3*x*sin(x);
    double fxU= (a1+gsl_vector_get(dc,0)) + (a2+gsl_vector_get(dc,1))*x + (a3+gsl_vector_get(dc,2))*x*sin(x);
    double fxL= (a1-gsl_vector_get(dc,0)) + (a2-gsl_vector_get(dc,1))*x + (a3-gsl_vector_get(dc,2))*x*sin(x);
    printf("%lg\t%lg\t%lg\t%lg\n", x, fx, fxU, fxL);
  } printf("\n\n");


  /* Exercise C */
  /* n = 10, already defined */
  /* Define x,y,dy like in exercise a */
  for(int i=0; i<n; i++) {
    gsl_vector_set(x,  i, xdata[i]);
    gsl_vector_set(y,  i, ydata[i]);
    gsl_vector_set(dy, i, dydata[i]);
  }

  /* m = 3, already defined */
  /* gsl_vector* c  = gsl_vector_alloc(m); */
  /* gsl_vector* dc = gsl_vector_alloc(m); */
  /* gsl_matrix* S  = gsl_matrix_alloc(m,m); */

  lssvd(x,y,dy,m,funcA,c,S);

  fprintf(stderr, "Exercise C fitting coefficients\n");
  for(int i=0; i<m; i++) gsl_vector_set(dc, i, sqrt(gsl_matrix_get(S,i,i)));
  
  /* Print result */
  fprintf(stderr,"f(x) = a1/x + a2 + a3*x\n[a1,a2,a3] = \n");
  for(int i=0; i<m; i++) fprintf(stderr,"%7.5lg +/- %7.5lg\n", gsl_vector_get(c,i), gsl_vector_get(dc,i)); 
  fprintf(stderr,"\n");
 
  /* Print function for plotting */
  a1 = gsl_vector_get(c,0), a2 = gsl_vector_get(c,1), a3 = gsl_vector_get(c,2);
  /* double x0=0.07, x1=3.0, dx=0.05; */
  for(double x=x0; x<x1; x+=dx) {

    double fx = a1/x  + a2  + a3*x;
    double fxU= (a1+gsl_vector_get(dc,0))/x + (a2+gsl_vector_get(dc,1)) + (a3+gsl_vector_get(dc,2))*x;
    double fxL= (a1-gsl_vector_get(dc,0))/x + (a2-gsl_vector_get(dc,1)) + (a3-gsl_vector_get(dc,2))*x;
    printf("%lg\t%lg\t%lg\t%lg\n", x, fx, fxU, fxL);
  } printf("\n\n");

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(dy);
  gsl_vector_free(c);
  gsl_vector_free(dc);
}
