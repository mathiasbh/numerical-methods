
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>
#include "header.h"

/* n by m (n>m) matrix A */
/* n by m orthogonal matrix U */
/* m by m square matrix S diagonal, non-negative real */
/* m by m orthogonal matrix V */
/* SVD: solve linear problem Ac=b, USV^Tc = b */


void svd(gsl_matrix* A, gsl_matrix* U, gsl_matrix* S, gsl_matrix* V, gsl_matrix* Sigma) {
  /* SVD of tall matrix, A^T A = VDV^T
     V: Eigenvectors
     D: diagonal matrix with eigenvalues */

 /* Input: A
    Output: S, U Sigma, V */
  assert(A->size2 == S->size2);
  assert(A->size2 == V->size2);
  assert(S->size1 == S->size2);
  assert(V->size1 == V->size2);
  assert(A->size1 == U->size1);
  assert(A->size2 == U->size2);

  int n = A->size1, m = A->size2;

  gsl_vector* e = gsl_vector_alloc(m);
  gsl_matrix* ATA = gsl_matrix_alloc(m,m);
  /* A^T * A */
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, ATA);
  

  jacobi_sweep(ATA,e,V);
  
  /* Construct S = D^(1/2) and D(-1/2) */
  gsl_matrix* Dinvsqrt = gsl_matrix_alloc(m,m);
  for(int i=0; i<m; i++) {
    gsl_matrix_set(S, i, i, sqrt(gsl_vector_get(e,i))); /* D^(1/2) */
    gsl_matrix_set(Dinvsqrt, i, i, 1.0/sqrt(gsl_vector_get(e,i))); /* D^(-1/2) */
  }

  /* Construct U = A V D^(-1/2) */
  gsl_matrix* AV = gsl_matrix_alloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, V, 0.0, AV);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AV, Dinvsqrt, 0.0, U);

  /* Construct covariance matrix Sigma = V S^-2 V^T = V D V^T */
  gsl_matrix* Dinv = gsl_matrix_alloc(m,m);
  gsl_matrix* VD = gsl_matrix_alloc(m,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Dinvsqrt, Dinvsqrt, 0.0, Dinv);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, Dinv, 0.0, VD); 
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, VD, V, 0.0, Sigma); 
  

  gsl_vector_free(e);
  gsl_matrix_free(ATA);
  gsl_matrix_free(Dinvsqrt);
  gsl_matrix_free(AV);
  gsl_matrix_free(Dinv);
  gsl_matrix_free(VD);
}

void svd_solve(gsl_matrix* U, gsl_matrix* S, gsl_matrix* V, gsl_vector* b, gsl_vector* c) {
  assert(S->size1 == S->size2);
  assert(V->size1 == V->size2);
  assert(V->size1 == U->size2);
  assert(S->size2 == U->size2);
  
  /* Back substitution */
  /* Solve Ax=b, or in this case SV^T*c = Sy = U^T b */
  gsl_vector* UTb = gsl_vector_alloc(U->size2);
  gsl_blas_dgemv(CblasTrans, 1.0, U, b, 0.0, UTb); /* U^T * b */

  int n = S->size2;
  for(int i=n-1; i>=0; i--) {
    double SumSy = 0;
    for(int k=i+1; k<n; k++) {
      SumSy -= gsl_matrix_get(S, i, k) * gsl_vector_get(UTb, k);
    }
    gsl_vector_set(UTb, i, (gsl_vector_get(UTb,i) + SumSy)/gsl_matrix_get(S, i, i)); /* y stored in UTb */
  }
  
  /* c = V*y */
  gsl_blas_dgemv(CblasNoTrans, 1.0, V, UTb, 0.0, c); /* U^T * b */

  gsl_vector_free(UTb);
}


void lssvd(gsl_vector* x, gsl_vector* y, gsl_vector* dy, 
	   int nf, double f(int k, double x), gsl_vector* c, gsl_matrix* Sigma) {
  assert(x->size == y->size);
  assert(x->size == dy->size);
 /* Input: x,y,dy, nf(number of parameters), function to fit, Sigma
    Output: c, solved from Ac = b */
  int n = x->size;
  int m = nf;

  gsl_matrix* A = gsl_matrix_alloc(n,m);
  gsl_vector* b = gsl_vector_alloc(n);
  gsl_matrix* V = gsl_matrix_alloc(m,m);
  gsl_matrix* S = gsl_matrix_alloc(m,m);
  gsl_matrix* U = gsl_matrix_alloc(n,m);

  /* Calculate bi = yi/dyi and Aik = fk(xi)/dyi */
  for(int i=0; i<n; i++) {
    gsl_vector_set(b, i, gsl_vector_get(y,i)/gsl_vector_get(dy,i));
    for(int k=0; k<m; k++) {
      gsl_matrix_set(A, i, k, f(k, gsl_vector_get(x,i))/gsl_vector_get(dy,i));
    }
  }

  /* Decomposition using singular value */
  svd(A,U,S,V,Sigma); /* Returns constructed U,S,V,Sigma */
  svd_solve(U,S,V,b,c);
  

  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_matrix_free(V);
  gsl_matrix_free(U);
  gsl_matrix_free(S);
}
