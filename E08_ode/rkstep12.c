/* Runge-Kutta stepper one two */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void rkstep12(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* dy) {
  int n = y0->size;
  gsl_vector* k0  = gsl_vector_alloc(n);
  gsl_vector* yt  = gsl_vector_alloc(n);
  gsl_vector* k12 = gsl_vector_alloc(n);
  f(t,y0,k0);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i) + gsl_vector_get(k0,i)*h/2.0);
  }
  
  f(t+h/2.0,yt,k12);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yh, i, gsl_vector_get(y0,i) + gsl_vector_get(k12,i)*h);
  }
  for(int i=0; i<n; i++) {
    gsl_vector_set(dy, i, (gsl_vector_get(k0,i) - gsl_vector_get(k12,i))*h/2.0);
  }
  gsl_vector_free(k0);
  gsl_vector_free(yt);
  gsl_vector_free(k12);
}
