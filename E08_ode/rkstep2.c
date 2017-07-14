/* Runge-Kutta stepper 2 - not embedded*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void rkstep2(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* dy) {
  int n = y0->size;
  int p = 2;
  gsl_vector* k1  = gsl_vector_alloc(n);
  gsl_vector* k2  = gsl_vector_alloc(n);
  gsl_vector* yt  = gsl_vector_alloc(n);
  
  f(t,y0,k1);
  for(int i=0; i<n; i++) gsl_vector_set(yt, i, gsl_vector_get(y0,i) + gsl_vector_get(k1,i)*h/2.0);
  
  f(t+h/2.0,yt,k2);
  for(int i=0; i<n; i++) gsl_vector_set(yh, i, gsl_vector_get(y0,i) + gsl_vector_get(k2,i)*h);

  for(int i=0; i<n; i++) {
    gsl_vector_set(dy, i, (gsl_vector_get(k1,i)*h
			   - gsl_vector_get(k1,i)*h/2.0
			   - gsl_vector_get(k2,i)*h/2.0)/(pow(2.0,p)-1.0));
  }
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(yt);
}
