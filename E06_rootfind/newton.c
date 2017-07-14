/* 
   Newton's method with back-tracking linesearch -
   derivatives in the Jacobian matrix are evaluated numerically.

   function f(x,fx): takes the input vector x, calculates the vector f(x) and puts it into vector fx;
   vector xstart: the starting point, upon exit it becomes the latest approximation to the root;
   double dx: the delta x to be used in numerical evaluation of the Jacobian matrix;
   double epsilon: the accuracy goal.
*/

#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "qrgivens.h"

void newton_rootfind(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double epsilon) {
  int n = x->size;
  gsl_matrix* jacobi = gsl_matrix_alloc(n,n); 
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* fx = gsl_vector_alloc(n);
  gsl_vector* z = gsl_vector_alloc(n);
  gsl_vector* fz = gsl_vector_alloc(n);
  int numberSteps = 0;


  do {
    numberSteps++;
    f(x,fx); // Calculate f(x) and store in fx
    for(int j=0; j<n; j++) {      /* eqn 7: numerical eval of jacobian */
      gsl_vector_set(x, j, gsl_vector_get(x, j) + dx);
      f(x,fz);
      gsl_vector_sub(fz,fx); // fz-fx
      for(int i=0; i<n; i++) gsl_matrix_set(jacobi, i, j, gsl_vector_get(fz, i)/dx);
      gsl_vector_set(x, j, gsl_vector_get(x, j) - dx); 
    }
    /* Solve J Dx = -f(x) */
    qr_givens_decomp(jacobi);
    gsl_vector_scale(fx, -1);
    qr_givens_solve(jacobi, fx, deltax); // solve J*deltax = -f(x) for deltax, eqn 5

    double lambda = 2.0; 
    double normfx = gsl_blas_dnrm2(fx);

    do { /* back-tracking line search */
      lambda/=2.0;
      gsl_vector_memcpy(z,x); // z <- x + lambda*deltax
      gsl_blas_daxpy(lambda, deltax, z); 
      f(z,fz);
    } while(gsl_blas_dnrm2(fz) > (1.0 - lambda/2.0)*normfx && lambda > 1.0/128.0);
    /* Copy z into x and fz into fx */
    gsl_vector_memcpy(x,z); gsl_vector_memcpy(fx,fz);
  } while(gsl_blas_dnrm2(fx) > epsilon && gsl_blas_dnrm2(deltax) > dx);

  printf("Newton root find steps: %i\n", numberSteps);

  gsl_matrix_free(jacobi);
  gsl_vector_free(z);
  gsl_vector_free(deltax);
  gsl_vector_free(fx);
  gsl_vector_free(fz);
}


/* User specifies jacobian */
void newton_rootfind_jacobi(void f(gsl_vector* x, gsl_vector* fx), 
			    gsl_matrix* jacobian(gsl_vector* x), gsl_vector* x, double dx, double epsilon) {
  int n = x->size;
  gsl_matrix* jacobi = gsl_matrix_alloc(n,n); 
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* fx = gsl_vector_alloc(n);
  gsl_vector* z = gsl_vector_alloc(n);
  gsl_vector* fz = gsl_vector_alloc(n);
  int numberSteps = 0;


  do {
    numberSteps++;
    f(x,fx); // Calculate f(x) and store in fx
    jacobi = jacobian(x);
    /* Solve J Dx = -f(x) */
    qr_givens_decomp(jacobi);
    gsl_vector_scale(fx, -1);
    qr_givens_solve(jacobi, fx, deltax); // solve J*deltax = -f(x) for deltax, eqn 5

    double lambda = 2.0; 
    double normfx = gsl_blas_dnrm2(fx);

    do { /* back-tracking line search */
      lambda/=2.0;
      gsl_vector_memcpy(z,x); // z <- x + lambda*deltax
      gsl_blas_daxpy(lambda, deltax, z); 
      f(z,fz);
    } while(gsl_blas_dnrm2(fz) > (1.0 - lambda/2.0)*normfx && lambda > 1.0/128.0);
    /* Copy z into x and fz into fx */
    gsl_vector_memcpy(x,z); gsl_vector_memcpy(fx,fz);
  } while(gsl_blas_dnrm2(fx) > epsilon && gsl_blas_dnrm2(deltax) > dx);

  printf("Newton root find steps: %i\n", numberSteps);

  gsl_matrix_free(jacobi);
  gsl_vector_free(z);
  gsl_vector_free(deltax);
  gsl_vector_free(fx);
  gsl_vector_free(fz);
}
