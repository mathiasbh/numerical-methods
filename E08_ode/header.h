#include <gsl/gsl_vector.h>

void driver(double* t, double b, double* h, gsl_vector* y, double acc, double eps,
	void stepper(double t, double h, gsl_vector *y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* err),
	void f(double t,gsl_vector*y,gsl_vector*dydt));


int driver_list(double* tlist, double b, double* h, gsl_vector* y, double** ylist, double acc, double eps, int max,
		 void stepper(double t, double h, gsl_vector *y0,
			      void f(double t, gsl_vector* y, gsl_vector* dydt),
			      gsl_vector* yh, gsl_vector* dy),
		void f(double t,gsl_vector*y,gsl_vector*dydt));


void rkstep12(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* err);

void rkstep45(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	      gsl_vector* yh, gsl_vector* dy);

void rkstep2(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	     gsl_vector* yh, gsl_vector* dy);

void rkstep5(double t, double h, gsl_vector* y0,
	      void f(double t, gsl_vector* y, gsl_vector* dydt),
	     gsl_vector* yh, gsl_vector* dy);
