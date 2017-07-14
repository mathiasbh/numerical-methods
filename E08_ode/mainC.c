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


int main() {
  int n = 2, max = 1000;
  double h=0.1, acc=1e-3, eps=1e-3, t=0, b=6.28;
  
  double* tlist = (double* )calloc(max,sizeof(double));
  double**ylist = (double**)calloc(max,sizeof(double*));
  gsl_vector* y = gsl_vector_alloc(n);

  printf("#RK12\n");
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  tlist[0]=t;
  
  int k12 = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep12,sinefunction);
  fprintf(stderr, "RK12: steps=%i\n", k12);
  for(int i=0; i<k12; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");

  printf("#RK45\n");
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  tlist[0]=t;
  
  int k45 = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep45,sinefunction);
  fprintf(stderr, "RK45: steps=%i\n", k45);
  for(int i=0; i<k45; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");

  printf("#RK2\n");
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  tlist[0]=t;
  
  int k2 = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep2,sinefunction);
  fprintf(stderr, "RK2 : steps=%i\n", k2);
  for(int i=0; i<k2; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");

  
  printf("#RK5\n");
  gsl_vector_set(y, 0, 0);
  gsl_vector_set(y, 1, 1);
  for(int i=0; i<max; i++) ylist[i] = (double*)calloc(n,sizeof(double));
  tlist[0]=t;
  
  int k5 = driver_list(tlist,b,&h,y,ylist,acc,eps,max,rkstep5,sinefunction);
  fprintf(stderr, "RK5 : steps=%i\n", k5);
  for(int i=0; i<k5; i++) printf("%lg %lg %lg\n", tlist[i], ylist[i][0], ylist[i][1]);
  printf("\n\n");
  
  gsl_vector_free(y);
  return 0;
}
