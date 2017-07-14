/* Some text */


#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>

void   qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
void   qr_gs_solve (const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x);
double qr_gs_absdet(const gsl_matrix* R);
void   qr_gs_invert(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* Ainverse);

void printm(const gsl_matrix* A);
void printv(const gsl_vector* v);
void matrix_set_rand(gsl_matrix* A);
void vector_set_rand(gsl_vector *v);


int main() {


  /* Allocate memory to matrices and vectors */
  int n = 4; /* Row */
  int m = 3; /* Column */
  gsl_matrix* A1= gsl_matrix_calloc(n,m); /* 1: Gram Schmidt non-square */
  gsl_matrix* R1= gsl_matrix_calloc(m,m);
  gsl_vector* b1= gsl_vector_calloc(m);
  gsl_vector* x1= gsl_vector_calloc(m);


  /* Set random matrix A non-square */
  matrix_set_rand(A1);

  /* Set random vector b */
  vector_set_rand(b1);

  /* ############################################################ */
  /* Test results - Gram Schmidt */
  printf("### Gram Schmidt test ###\nInput system (non-square matrix)\n\n");
  printf("Matrix A:\n"); printm(A1);
  printf("Vector b:\n"); printv(b1);

  qr_gs_decomp(A1, R1); /* turn A into Q */
  printf("Matrix Q:\n"); printm(A1);
  printf("Matrix R:\n"); printm(R1);

  gsl_matrix* QTQ = gsl_matrix_calloc(m,m);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A1, A1, 0.0, QTQ);
  printf("Matrix Q^T*Q = I\n"); printm(QTQ);

  gsl_matrix* QR = gsl_matrix_calloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A1, R1, 0.0, QR);
  printf("Matrix Q*R - should be A:\n"); printm(QR);

  /* Square matrix */
  /* Timing */
  clock_t start_t, end_t;
  start_t = clock();
  n = 10; m = 10;
  gsl_matrix* A2 = gsl_matrix_calloc(n,m); /* 2: Gram Schmidt square/linear system */
  gsl_vector* b2 = gsl_vector_calloc(m);
  gsl_matrix* R2 = gsl_matrix_calloc(m,m);
  gsl_vector* x2 = gsl_vector_calloc(m);

  /* Set random matrix A - square */
  matrix_set_rand(A2); vector_set_rand(b2);

  printf("Linear system\nInput system (square matrix)\n\n");
  printf("Matrix A:\n"); printm(A2);
  printf("Vector b:\n"); printv(b2);

  qr_gs_decomp(A2, R2);
  qr_gs_solve(A2, R2, b2, x2);
  printf("Result, vector x:\n"); printv(x2);
  
  /* Form A = Q*R for further test */
  gsl_matrix* QR2 = gsl_matrix_calloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A2, R2, 0.0, QR2);

  gsl_vector* b_test = gsl_vector_calloc(n);
  gsl_blas_dgemv(CblasNoTrans, 1.0, QR2, x2, 0.0, b_test);
  printf("Vector A*x - should be b:\n"); printv(b_test);


  printf("|det(A)| = %7.3lg\n", qr_gs_absdet(R2));

  gsl_matrix* Ainverse = gsl_matrix_alloc(n,m);
  qr_gs_invert(A2, R2, Ainverse);
  printf("Matrix A^-1:\n"); printm(Ainverse);
  gsl_matrix* AinvA = gsl_matrix_alloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ainverse, QR2, 0.0, AinvA);
  printf("Matrix A^-1*A = I:\n"); printm(AinvA);
  
  /* Free everything */
  gsl_matrix_free(A1);
  gsl_matrix_free(R1);
  gsl_matrix_free(AinvA);
  gsl_matrix_free(Ainverse);
  gsl_matrix_free(QR2);
  gsl_matrix_free(A2);
  gsl_matrix_free(R2);
  gsl_matrix_free(QTQ);
  gsl_matrix_free(QR);

  gsl_vector_free(b1);
  gsl_vector_free(x1);
  gsl_vector_free(b_test);
  gsl_vector_free(b2);
  gsl_vector_free(x2);


  end_t = clock();
  double diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> Gram Schmidt timing: %f ms <<<<\n\n", diff_t*1000);
  return 0;
}

