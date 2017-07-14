
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "header.h"

void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, 
	   int nf, double f(int k, double x), gsl_vector* c, gsl_matrix* S) {
  int n = x->size;
  int m = nf;
  
  gsl_matrix* A = gsl_matrix_alloc(n,m);
  gsl_matrix* R = gsl_matrix_alloc(m,m);
  gsl_matrix* Rinv = gsl_matrix_alloc(m,m);
  gsl_vector* b = gsl_vector_alloc(n);
  gsl_vector* p = gsl_vector_alloc(n);

  /* Calculate bi = yi/dyi and Aik = fk(xi)/dyi */
  for(int i=0; i<n; i++) {
    gsl_vector_set(b, i, gsl_vector_get(y,i)/gsl_vector_get(dy,i));
    for(int k=0; k<m; k++) {
      gsl_matrix_set(A, i, k, f(k, gsl_vector_get(x,i))/gsl_vector_get(dy,i));
    }
  }
  
  /* Decomposition using Givens rotation */
  qr_givens_decomp(A);
  qr_givens_unpackR(A, R);
  qr_givens_solve(A, b, p);
  for (int i=0; i<m; i++) {
    gsl_vector_set(c, i, gsl_vector_get(p, i));
  }

  /* Calculate estimates of errors Delta c
     Delta c = sqrt( (R^T * R)^-1 )|diagonal elements */
  qr_givens_invert(R, Rinv);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Rinv, Rinv, 0.0, S);

  
  
  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(Rinv);
  gsl_vector_free(b);
  gsl_vector_free(p);
}
