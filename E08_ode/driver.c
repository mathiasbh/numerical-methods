#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void driver(double* t, double b, double* h, gsl_vector* y, double acc, double eps,
	void stepper(double t, double h, gsl_vector *y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* dy),
	void f(double t,gsl_vector*y,gsl_vector*dydt)) {
  int n = y->size;
  int k = 0;
  gsl_vector* yh = gsl_vector_alloc(n);
  gsl_vector* dy = gsl_vector_alloc(n);
  double err, normy, tol;
  double a = *t;
  
  if(*h > (b-a)) *h = (b-a);

  while(*t < b) {
    stepper(*t,*h,y,f,yh,dy);
    err = gsl_blas_dnrm2(dy);
    normy = gsl_blas_dnrm2(yh);
    tol = (normy*eps+acc)*sqrt(*h/(b-a));
    if(err<tol) { /* Accept step */
      k++;
      *t+= *h;
      gsl_vector_memcpy(y, yh);
    }
    if(err>0) *h*=pow(tol/err,0.25)*0.95;
    else *h*=2;
  }
  gsl_vector_free(yh);
  gsl_vector_free(dy);
}

int driver_list(double* tlist, double b, double* h, gsl_vector* y, double** ylist, double acc, double eps, int max,
		 void stepper(double t, double h, gsl_vector *y0,
		     void f(double t, gsl_vector* y, gsl_vector* dydt),
		     gsl_vector* yh, gsl_vector* dy),
		 void f(double t,gsl_vector*y,gsl_vector*dydt)) {
  
  int n = y->size;
  int k = 0;
  gsl_vector* yh = gsl_vector_alloc(n);
  gsl_vector* dy = gsl_vector_alloc(n);
  double err, normy, tol, t;
  double a = tlist[0];
  
  if(*h > (b-a)) *h = (b-a);

  while(tlist[k] < b) {
    t = tlist[k];
    stepper(t,*h,y,f,yh,dy);
    err = gsl_blas_dnrm2(dy);
    normy = gsl_blas_dnrm2(yh);
    tol = (normy*eps+acc)*sqrt(*h/(b-a));
    if(err<tol) { /* Accept step */
      k++;
      if(k>max-1) return k;
      tlist[k] = t + *h;
      gsl_vector_memcpy(y, yh);
      for(int i=0; i<n; i++) ylist[k][i] = gsl_vector_get(y,i);
    }
    if(err>0) *h*=pow(tol/err,0.25)*0.95;
    else *h*=2;
  }
  
  gsl_vector_free(yh);
  gsl_vector_free(dy);
  return k;
}
