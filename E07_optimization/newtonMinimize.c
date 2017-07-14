#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "qrgivens.h"

void newton_minimize(void f(gsl_vector* x, double* fx),
		     void grad(gsl_vector* x, gsl_vector* df),
		     void hess(gsl_vector* x, gsl_matrix* H),
		     gsl_vector* x, double dx, double epsilon) {
  int n = x->size;
  gsl_matrix* H = gsl_matrix_alloc(n,n);
  gsl_vector* z = gsl_vector_alloc(n);
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* df = gsl_vector_alloc(n);
  double fx, fz;
  int numberSteps = 0;

  do {
    numberSteps++;
    f(x, &fx); // Calculate f(x) and store in fx, same for gradientf and Hessian
    grad(x, df);
    hess(x, H);
    
    qr_givens_decomp(H);
    gsl_vector_scale(df, -1);
    qr_givens_solve(H, df, deltax); // solve H*deltax = -df(x) for deltax, eqn 5

    double lambda = 2.0, alpha = 1e-4, condition; // condition: condition for while loop

    do { /* back-tracking line search */
      lambda/=2.0;
      gsl_vector_memcpy(z,x); // z <- x + lambda*deltax
      gsl_blas_daxpy(lambda, deltax, z); 
      f(z,&fz);
      
      /* Equation 8: condition */
      gsl_blas_ddot(deltax, df, &condition);
      condition *= alpha*lambda;
      condition += fx;
    } while(fabs(fz) > condition && lambda > 1.0/128.0);
    
    /* Copy z into x */
    gsl_vector_memcpy(x,z); fx = fz;
  } while(gsl_blas_dnrm2(deltax) > dx && gsl_blas_dnrm2(df) > epsilon);

  printf("Newton minimization steps: %i\n", numberSteps);

  gsl_matrix_free(H);
  gsl_vector_free(z);
  gsl_vector_free(deltax);
  gsl_vector_free(df);

}


/* Quasi newton method with Broyden's update */
void newton_minimize_quasi(void f(gsl_vector* x, double* fx),
		     void grad(gsl_vector* x, gsl_vector* df),
		     gsl_vector* x, double dx, double epsilon) {
  int n = x->size;
  gsl_matrix* Hinv = gsl_matrix_alloc(n,n);
  gsl_matrix* dH = gsl_matrix_alloc(n,n); // update Hessian
  gsl_matrix* Hplaceholder = gsl_matrix_alloc(n,n);
  gsl_vector* z = gsl_vector_alloc(n);
  gsl_vector* deltax = gsl_vector_alloc(n);
  gsl_vector* df = gsl_vector_alloc(n);
  gsl_vector* dfz = gsl_vector_alloc(n);
  gsl_vector* s = gsl_vector_alloc(n);
  gsl_vector* y = gsl_vector_alloc(n);
  gsl_vector* u = gsl_vector_alloc(n);
  double denom = 0;
  double fx, fz;
  int numberSteps = 0;

  f(x,&fx);
  grad(x, df);
  
  gsl_matrix_set_identity(Hinv); // First iteration of inverse Hessian


  /* f(x, &fx); // Calculate f(x) and store in fx, same for gradientf and Hessian */
  
  do {
    numberSteps++;

    /* deltax = -H^- * grad fx */
    gsl_blas_dgemv(CblasNoTrans, -1.0, Hinv, df, 0.0, deltax);
    double lambda = 2.0, alpha = 1e-4, condition; // condition for while loop
    gsl_vector_memcpy(s,deltax); gsl_vector_scale(s,lambda);

    do { /* back-tracking line search */
      lambda/=2.0; gsl_vector_scale(s,lambda);
      gsl_vector_memcpy(z,x); // z <- x + lambda*deltax
      gsl_blas_daxpy(lambda, deltax, z);
      f(z,&fz);
      
      /* Equation 8: condition */
      gsl_blas_ddot(deltax, df, &condition);
      condition *= alpha*lambda;
      condition += fabs(fx); // f(x+s) < f(x) + alpha *lambda* Dx^T * df(x)
      if(fabs(fz) < condition) break;
      if(gsl_blas_dnrm2(s) < dx) {gsl_matrix_set_identity(Hinv); break;}
    } while(1); //while(fabs(fz) > condition && lambda > 1.0/512.0);
    
    /* Update H^-1 */
    /* y = df(x+s) - df(x) */
    grad(z,dfz); gsl_vector_memcpy(y,dfz); gsl_blas_daxpy(-1.0, df, y);

    /* Calculate update, store in s and dH */
    /* u -> s - H^-1 * y, where s = lambda*deltax */
    gsl_vector_memcpy(u,s);
    gsl_blas_dgemv(CblasNoTrans, -1.0, Hinv, y, 1.0, u);
    
    /* dH -> (s-H^-1*y) * s^T = u*s^T, vector*vector^T = matrix, no blas function for this */
    for(int i=0; i<n; i++) for(int j=0; j<n; j++) {
	gsl_matrix_set(dH,i,j, gsl_vector_get(u,i)*gsl_vector_get(s,j));
      }

    /* (s-H^-1*y) * s^T * H^-1 = dH*H^-1, matrix*matrix: numerator stored in Hplaceholder */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, dH, Hinv, 0.0, Hplaceholder);

    /* H^-1 * s, use u as placeholder */
    gsl_blas_dgemv(CblasNoTrans, 1.0, Hinv, s, 0.0, u);

    /* y^T H^-1 * s = y^T*u, stored in denom */
    gsl_blas_ddot(y, u, &denom);
    /* fprintf(stderr, "denom = %lg\n", denom); */
    /* Calculate update: H^-1 = H^-1 * (s - H^-1*y) s^T H^-1 / y^T H^-1 s */
    gsl_matrix_scale(Hplaceholder, 1.0/denom);
    gsl_matrix_add(Hinv, Hplaceholder);

    /* Copy z into x */
    gsl_vector_memcpy(x,z); fx = fz; gsl_vector_memcpy(df,dfz);
  } while(gsl_blas_dnrm2(deltax) > dx && gsl_blas_dnrm2(df) > epsilon);

  printf("Newton quasi minimization steps: %i\n", numberSteps);

  gsl_matrix_free(Hinv);
  gsl_matrix_free(dH);
  gsl_matrix_free(Hplaceholder);
  gsl_vector_free(z);
  gsl_vector_free(deltax);
  gsl_vector_free(df);
  gsl_vector_free(dfz);
  gsl_vector_free(s);
  gsl_vector_free(y);
  gsl_vector_free(u);
}
