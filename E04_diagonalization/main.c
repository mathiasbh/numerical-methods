
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "jacobi.h"
#include <time.h>

void printing(gsl_matrix* A, gsl_vector* e, gsl_matrix* V) {
  printm(A);
  printf("Eigenvalues:\n"); printv(e);
  printf("Eigenvectors:\n"); printm(V);
}

int main() {
  /* Form test matrix */
  int n = 6;
  gsl_matrix* A = gsl_matrix_alloc(n,n);
  gsl_matrix* V = gsl_matrix_alloc(n,n);
  gsl_vector* e = gsl_vector_alloc(n);

  /* Random function modified to produce symmetric matrix */
  matrix_set_rand(A);
  printf("Matrix A (real symmetric) used for testing:\n\n\n"); printm(A);

  
  clock_t tic = clock();
  printf("Test cyclic sweep diagonalization\n");
  
  /* Run jacobi diagonalization using sweep */
  int rot1 = jacobi_sweep(A, e, V);
  printf("Matrix A (Jacobi'ed):\n"); 
  printing(A,e,V);

  clock_t toc = clock();
  printf("Number of rotations: %d\n", rot1);
  printf("Elapsed time: %f ms\n\n\n", (double)(toc-tic)*1000/CLOCKS_PER_SEC);

  /* printf("Restore A test\n"); */
  restore_matrix(A); 
  /* printm(A); */

  /* Run jacobi diagonalization using sweep row by row */
  tic = clock();

  printf("Test row by row diagonalization\n");
  int rot2 = jacobi_sweep_row(A, e, V, n-1, 1);
  printf("Matrix A (Jacobi'ed):\n"); 
  printing(A,e,V);

  toc = clock();
  printf("Number of rotations: %d\n", rot2);
  printf("Elapsed time: %f ms\n\n\n", (double)(toc-tic)*1000/CLOCKS_PER_SEC);


  restore_matrix(A);
  /* Run jacobi diagonalization using sweep row by row eliminating largest first */
  tic = clock();

  printf("Test row by row - eliminate largest first - diagonalization\n");
  int rot3 = jacobi_sweep_row_max(A, e, V, n-1, 1);
  printf("Matrix A (Jacobi'ed):\n"); 
  printing(A,e,V);

  toc = clock();
  printf("Number of rotations: %d\n", rot3);
  printf("Elapsed time: %f ms\n\n\n", (double)(toc-tic)*1000/CLOCKS_PER_SEC);
  

  /* Free memory */
  gsl_matrix_free(A); gsl_matrix_free(V); gsl_vector_free(e);
  return 0;
}
