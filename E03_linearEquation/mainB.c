
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>

void   qr_givens_decomp(gsl_matrix* A);
void   qr_givens_solve (gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
double qr_givens_absdet(const gsl_matrix* QR);
void   qr_givens_invert(const gsl_matrix* QR, gsl_matrix* Ainverse); 
void   qr_givens_unpackQ(gsl_matrix* QR, gsl_matrix *Q);
void   qr_givens_unpackR(gsl_matrix* QR, gsl_matrix *R);


void printm(const gsl_matrix* A);
void printv(const gsl_vector* v);
void matrix_set_rand(gsl_matrix* A);
void vector_set_rand(gsl_vector *v);

int main() {

  /* Allocate memory to matrices and vectors */
  int n = 4; /* Row */
  int m = 3; /* Column */
  gsl_matrix* A1= gsl_matrix_calloc(n,m);
  gsl_matrix* R1= gsl_matrix_calloc(m,m);
  gsl_vector* b1= gsl_vector_calloc(m);
  gsl_vector* x1= gsl_vector_calloc(m);
  gsl_matrix* Q1= gsl_matrix_alloc(n,m);
  gsl_matrix* A1_remade = gsl_matrix_calloc(n,m);

  matrix_set_rand(A1); vector_set_rand(b1);
  /* ############################################################ */
  /* Test results - Givens */
  printf("### Givens rotation test ###\nInput system\n\n");


  printf("Matrix A:\n"); printm(A1);
  printf("Vector b:\n"); printv(b1);

  qr_givens_decomp(A1); /* Turn A into QR/ modified R */
  printf("Matrix decomposition (R with theta stored lower matrix):\n"); printm(A1);
  
  /* Unpack Q and R */
  qr_givens_unpackQ(A1, Q1);
  qr_givens_unpackR(A1, R1);
  printf("Matrix Q unpacked:\n"); printm(Q1);
  printf("Matrix R unpacked:\n"); printm(R1);

  /* Q^T * Q */
  gsl_matrix* QTQ = gsl_matrix_calloc(m,m);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Q1, Q1, 0.0, QTQ);
  printf("Matrix Q^T*Q = I\n"); printm(QTQ);

  /* Form A = Q*R  */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Q1, R1, 0.0, A1_remade);
  printf("Matrix Q*R, should be A:\n"); printm(A1_remade);

  /* Square matrix */
  /* Timing */
  clock_t start_t, end_t;
  start_t = clock();
  printf("Linear system\nInput system (square matrix)\n\n");
  n = 10; m = 10;
  gsl_matrix* A2 = gsl_matrix_calloc(n,m); /* 2: Givens square/linear system */
  gsl_matrix* AA2= gsl_matrix_calloc(n,m);
  gsl_vector* b2 = gsl_vector_calloc(m);
  gsl_matrix* R2 = gsl_matrix_calloc(m,m);
  gsl_vector* x2 = gsl_vector_calloc(m);

  /* Set random matrix and vector b */
  matrix_set_rand(A2); vector_set_rand(b2);
  gsl_matrix_memcpy(AA2, A2); /* copy A */
  printf("Matrix A:\n"); printm(A2);
  printf("Vector b:\n"); printv(b2);
  
  qr_givens_decomp(A2); /* A into QR (or modified R with angles in lower part of matrix) */
  qr_givens_solve(A2, b2, x2);
  printf("Result vector x:\n"); printv(x2);

  gsl_vector* b2_test = gsl_vector_calloc(n);
  gsl_blas_dgemv(CblasNoTrans, 1.0, AA2, x2, 0.0, b2_test);
  printf("Check result A*x, should be b:\n"); printv(b2_test);

  double det2 = qr_givens_absdet(A2);
  printf("Determinant: %lg\n\n", det2);

  gsl_matrix* Ainverse2 = gsl_matrix_calloc(n,m);
  qr_givens_invert(A2, Ainverse2);
  printf("Matrix A^-1\n"); printm(Ainverse2);

  gsl_matrix* AinvA2 = gsl_matrix_alloc(n,m);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Ainverse2, AA2, 0.0, AinvA2);
  printf("Matrix A^-1*A = I:\n"); printm(AinvA2);  

  /* Free everything */  
  gsl_matrix_free(A1);
  gsl_matrix_free(R1);
  gsl_vector_free(b1);
  gsl_vector_free(x1);
  gsl_matrix_free(Q1);
  gsl_matrix_free(QTQ);
  gsl_matrix_free(A1_remade);
  gsl_matrix_free(A2);
  gsl_matrix_free(AA2);
  gsl_matrix_free(R2);
  gsl_vector_free(b2);
  gsl_vector_free(x2);
  gsl_vector_free(b2_test);  
  gsl_matrix_free(Ainverse2);
  gsl_matrix_free(AinvA2);

  end_t = clock();
  double diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> Givens timing: %f ms <<<<\n\n", diff_t*1000);
  return 0;
}

