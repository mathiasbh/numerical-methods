#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "header.h"

void sinefunction(double t, gsl_vector* y, gsl_vector* dydt) {
  gsl_vector_set(dydt,0,  gsl_vector_get(y,1));
  gsl_vector_set(dydt,1, -gsl_vector_get(y,0));
}

void testf(double t, gsl_vector* y, gsl_vector* dydt) {
  /* y'' + 2*y' + y = 0 */
  gsl_vector_set(dydt,0,  gsl_vector_get(y,1));
  gsl_vector_set(dydt,1, -2*gsl_vector_get(y,1) - gsl_vector_get(y,0));
}

int main() {
  int n = 2, max = 1000;
  double* tlist = (double* )calloc(max,sizeof(double));
  double**ylist = (double**)calloc(max,sizeof(double*));
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  double h=0.1, acc=1e-3, eps=1e-3, t=0, b=6.28;
  tlist[0]=t;
  
  int kOut = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep12,sinefunction);
  for(int i=0; i<kOut; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");


  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  h=0.1; t=0;
  tlist[0]=t;

  kOut = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep12,testf);
  for(int i=0; i<kOut; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");
  
  gsl_vector_free(y);
  return 0;
}
