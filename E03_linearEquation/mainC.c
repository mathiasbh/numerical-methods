
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>

void qr_lu_decomp(gsl_matrix* A, gsl_matrix* L, gsl_matrix* U);
void qr_lu_solve(gsl_matrix* L, gsl_matrix* U, gsl_vector* b, gsl_vector* x);
void qr_lu_invert(gsl_matrix* L, gsl_matrix* U, gsl_matrix* Ainverse);

void printm(const gsl_matrix* A);
void printv(const gsl_vector* v);
void matrix_set_rand(gsl_matrix* A);
void vector_set_rand(gsl_vector *v);

int main() {
  /* Timing */
  clock_t start_t, end_t;
  start_t = clock();

  /* Allocate memory to matrices and vectors */
  int n = 10; /* Row */
  int m = 10; /* Column */
  gsl_matrix* A= gsl_matrix_calloc(n,m);
  gsl_vector* b= gsl_vector_calloc(m);
  gsl_vector* x= gsl_vector_calloc(m);
  gsl_matrix* A_remade = gsl_matrix_calloc(n,m);
  gsl_matrix* L = gsl_matrix_calloc(n,n);
  gsl_matrix* U = gsl_matrix_calloc(n,n);
  gsl_vector* b_test = gsl_vector_calloc(n);
  gsl_matrix* Ainverse = gsl_matrix_calloc(n,m);
  gsl_matrix* AinvA = gsl_matrix_alloc(n,m);

  matrix_set_rand(A); vector_set_rand(b);

  /* ############################################################ */
  /* Test results - LU deomposition */
  printf("### LU decomposition test ###\nInput system (square matrix)\n\n");

  printf("Matrix A:\n"); printm(A);
  printf("Vector b:\n"); printv(b);

  
  qr_lu_decomp(A, L, U); /* A into QR (or modified R with angles in lower part of matrix) */
  printf("Matrix L:\n"); printm(L);
  printf("Matrix U:\n"); printm(U);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, L, U, 0.0, A_remade);
  printf("Matrix L*U = A\n"); printm(A_remade);
  
  qr_lu_solve(L, U, b, x);
  printf("Result vector x:\n"); printv(x);

  
  gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, b_test);
  printf("Check result A*x, should be b:\n"); printv(b_test);


  qr_lu_invert(L, U, Ainverse);
  printf("Matrix A^-1\n"); printm(Ainverse);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ainverse, A, 0.0, AinvA);
  printf("Matrix A^-1*A = I:\n"); printm(AinvA);

  /* Free everything */  
  gsl_matrix_free(A);
  gsl_matrix_free(A_remade);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(b_test);
  gsl_matrix_free(L);
  gsl_matrix_free(U);
  gsl_matrix_free(Ainverse);
  gsl_matrix_free(AinvA);

  end_t = clock();
  double diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> LU timing: %f ms <<<<\n\n ", diff_t*1000);
  return 0;
}

