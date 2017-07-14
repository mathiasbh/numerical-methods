#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "header.h"

void sinefunction(double t, gsl_vector* y, gsl_vector* dydt) {
  /* y'' + y = 0 */
  gsl_vector_set(dydt,0,  gsl_vector_get(y,1));
  gsl_vector_set(dydt,1, -gsl_vector_get(y,0));
}

void testf(double t, gsl_vector* y, gsl_vector* dydt) {
  /* y'' + 2*y' + y = 0 */
  gsl_vector_set(dydt,0,  gsl_vector_get(y,1));
  gsl_vector_set(dydt,1, -2*gsl_vector_get(y,1) - gsl_vector_get(y,0));
}

int main() {
  int n = 2;
  double t=0, b=1.570;
  double h=0.01, acc=1e-5, eps=1e-5;
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);

  printf("Solution should be sin(x)\n");
  printf("Initial: t=%lg, h=%lg, y(a)=", t, h);
  for(int i=0; i<y->size; i++) printf("%lg, ", gsl_vector_get(y,i));
  printf("\n");
  
  driver(&t,b,&h,y,acc,eps,rkstep12,sinefunction);

  printf("Final  : t=%lg, h=%lg, y(b)=", t, h);
  for(int i=0; i<y->size; i++) printf("%lg, ", gsl_vector_get(y,i));
  printf("\n\n");



  t = 0; h = 0.01;
  gsl_vector_set(y, 0, 1.0);
  gsl_vector_set(y, 1, 3.0);

  printf("Solution should be y(t)=a*exp(-t) + b*t*exp(-t)\n");
  printf("Initial: t=%lg, h=%lg, y(a)=", t, h);
  for(int i=0; i<y->size; i++) printf("%lg, ", gsl_vector_get(y,i));
  printf("\n");
  
  driver(&t,b,&h,y,acc,eps,rkstep12,testf);

  printf("Final  : t=%lg, h=%lg, y(b)=", t, h);
  for(int i=0; i<y->size; i++) printf("%lg, ", gsl_vector_get(y,i));
  printf("\n\n");
  
  gsl_vector_free(y);
  return 0;
}
