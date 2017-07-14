/* Runge-Kutta stepper four five */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

void rkstep5(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* dy) {
  int n = y0->size;
  int p = 5;
  gsl_vector* k1  = gsl_vector_alloc(n);
  gsl_vector* k2  = gsl_vector_alloc(n);
  gsl_vector* k3  = gsl_vector_alloc(n);
  gsl_vector* k4  = gsl_vector_alloc(n);
  gsl_vector* k5  = gsl_vector_alloc(n);
  gsl_vector* k6  = gsl_vector_alloc(n);
  gsl_vector* yt  = gsl_vector_alloc(n);
  
  f(t,y0,k1);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i)
		   + gsl_vector_get(k1,i)*h/4.0);
  }
  
  f(t+h/4.0,yt,k2);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i)
		   + gsl_vector_get(k1,i)*h*3/32.0
		   + gsl_vector_get(k2,i)*h*9/32.0);
  }

  f(t+h*3/8.0,yt,k3);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i)
		   + gsl_vector_get(k1,i)*h*1932/2197.0
		   - gsl_vector_get(k2,i)*h*7200/2197.0
		   + gsl_vector_get(k3,i)*h*7296/2197.0);
  }

  f(t+h*12/13.0,yt,k4);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i)
		   + gsl_vector_get(k1,i)*h*439/216.0
		   - gsl_vector_get(k2,i)*h*8.0
		   + gsl_vector_get(k3,i)*h*3680/513.0
		   - gsl_vector_get(k4,i)*h*845/4104.0);
  }

  f(t+h,yt,k5);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yt, i, gsl_vector_get(y0,i)
		   - gsl_vector_get(k1,i)*h*8/27.0
		   + gsl_vector_get(k2,i)*h*2.0
		   - gsl_vector_get(k3,i)*h*3544/2565.0
		   + gsl_vector_get(k4,i)*h*1859/4104.0
		   - gsl_vector_get(k5,i)*h*11/40.0);
  }

  f(t+h/2.0,yt,k6);
  for(int i=0; i<n; i++) {
    gsl_vector_set(yh, i, gsl_vector_get(y0,i)
		   + gsl_vector_get(k1,i)*h*16/135.0
		   + gsl_vector_get(k3,i)*h*6656/12825.0
		   + gsl_vector_get(k4,i)*h*28561/56430.0
		   - gsl_vector_get(k5,i)*h*9/50.0
		   + gsl_vector_get(k6,i)*h*2/55.0);

    gsl_vector_set(dy, i, (gsl_vector_get(k1,i)*h
			 - gsl_vector_get(k1,i)*h/2.0
			   - gsl_vector_get(k6,i)*h/2.0)/(pow(2,p)-1));
  }
  gsl_vector_free(k1);
  gsl_vector_free(k2);
  gsl_vector_free(k3);
  gsl_vector_free(k4);
  gsl_vector_free(k5);
  gsl_vector_free(k6);

  gsl_vector_free(yt);
}
