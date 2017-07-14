#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "header.h"

int inv_iter(gsl_matrix* A, gsl_vector* b, double s, double eps, int maxIter, double *eigenValue) {
  assert(A->size1 == A->size2); assert(A->size1 == b->size);
  int n = A->size1, step = 0;
  double lam, err, eigOld, eigNew = s;

  gsl_vector_scale(b,1.0/gsl_blas_dnrm2(b)); // initial guess normalized
  gsl_matrix* M = gsl_matrix_alloc(n,n); /* temp matrix */
  gsl_vector* y = gsl_vector_alloc(n); /* Store linear solve result */
  gsl_vector* c = gsl_vector_alloc(n); /* temp vector */
  for(int i=0; i<n; i++) gsl_matrix_set(A,i,i,gsl_matrix_get(A,i,i) - s);

  do {
    eigOld = eigNew;

    /* Calculate y[i+1] by solving linear system */
    gsl_matrix_memcpy(M,A);
    qr_givens_decomp(M);
    qr_givens_solve(M, b, y);

    /* Iterate: b = y/norm(y) */
    gsl_vector_scale(y,1.0/gsl_blas_dnrm2(y));
    gsl_vector_memcpy(b,y);

    /* Rayleigh quotient: estimate eigenvalue and error */
    gsl_blas_dgemv(CblasNoTrans, 1.0, A, b, 0.0, c); // A*b = c
    gsl_blas_ddot(b, c, &lam); // b'*A*b
    eigNew = s + lam;
    err = fabs(eigNew - eigOld);
    step++;
  } while(step < maxIter && err > eps);

  gsl_vector_memcpy(b,y);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(M);

  *eigenValue = eigNew;
  return step;
}


int rayleigh(gsl_matrix* A, gsl_vector* b, double s, double eps, int maxIter, double *eigenValue) {
  assert(A->size1 == A->size2); assert(A->size1 == b->size);
  int n = A->size1, step = 0;
  double lam, err, eigOld, eigNew = s;

  gsl_vector_scale(b,1.0/gsl_blas_dnrm2(b)); // initial guess normalized
  gsl_matrix* M = gsl_matrix_alloc(n,n); /* temp matrix */
  gsl_matrix* W = gsl_matrix_alloc(n,n); /* temp matrix */
  gsl_vector* y = gsl_vector_alloc(n); /* Store linear solve result */
  gsl_vector* c = gsl_vector_alloc(n); /* temp vector */ 
  
  do {
    /* Copy A to W, and make (A - s*I) dynamic for each step */
    gsl_matrix_memcpy(W,A);
    for(int i=0; i<n; i++) gsl_matrix_set(W,i,i,gsl_matrix_get(W,i,i) - eigNew);
    eigOld = eigNew;

    /* Calculate y[i+1] by solving linear system */
    gsl_matrix_memcpy(M,W);
    qr_givens_decomp(M);
    qr_givens_solve(M, b, y);

    /* Iterate: b = y/norm(y) */
    gsl_vector_scale(y,1.0/gsl_blas_dnrm2(y));
    gsl_vector_memcpy(b,y);

    /* Rayleigh quotient: estimate eigenvalue and error */
    gsl_blas_dgemv(CblasNoTrans, 1.0, W, b, 0.0, c); // A*b = c
    gsl_blas_ddot(b, c, &lam); // b'*A*b
    eigNew += lam/3; // update eigNew and thus s, division to avoid jumping to next eigenvalue

    err = fabs(eigNew - eigOld);
    step++;
  } while(step < maxIter && err > eps);

  gsl_vector_memcpy(b,y);
  gsl_vector_free(y);
  gsl_vector_free(c);
  gsl_matrix_free(M);
  gsl_matrix_free(W);
  *eigenValue = eigNew;
  return step;
}
