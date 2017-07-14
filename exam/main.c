#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include "header.h"

int main() {

  /* Allocate memory to matrices and vectors */
  int n = 3, m = 3, maxIter = 1000, steps;
  double epsilon = 1e-6, eigenValue;
  gsl_matrix* A   = gsl_matrix_alloc(n,m);
  gsl_matrix* B   = gsl_matrix_alloc(n,m);
  gsl_vector* x   = gsl_vector_alloc(m);
  gsl_vector* y   = gsl_vector_alloc(m);
  gsl_matrix*tmpA = gsl_matrix_alloc(n,m);
  gsl_vector*tmpx = gsl_vector_alloc(m);
  gsl_matrix*tmpB = gsl_matrix_alloc(n,m);
  gsl_vector*tmpy = gsl_vector_alloc(m);

  /* Set random vectors */
  vector_set_rand(x); vector_set_rand(y);
  
  /* Produce matrix A with known eigenvalues and eigenvectors */
  gsl_matrix_set(A,0,0, 1.5); gsl_matrix_set(A,0,1, 1   ); gsl_matrix_set(A,0,2, -0.5 );
  gsl_matrix_set(A,1,0, 2.5); gsl_matrix_set(A,1,1, 0.75); gsl_matrix_set(A,1,2, -1.25);
  gsl_matrix_set(A,2,0, 1.5); gsl_matrix_set(A,2,1, 0.75); gsl_matrix_set(A,2,2, -0.25);
  gsl_matrix_memcpy(tmpA,A); gsl_vector_memcpy(tmpx,x);


  gsl_matrix_set(B,0,0, 2); gsl_matrix_set(B,0,1, 3); gsl_matrix_set(B,0,2, 1 );
  gsl_matrix_set(B,1,0, 0); gsl_matrix_set(B,1,1, 4); gsl_matrix_set(B,1,2, 0);
  gsl_matrix_set(B,2,0, 2); gsl_matrix_set(B,2,1, 0); gsl_matrix_set(B,2,2, 1);
  gsl_matrix_memcpy(tmpB,B); gsl_vector_memcpy(tmpy,y);
  
  /* Timing */
  clock_t start_t, end_t;
  start_t = clock();

  printf("Testing inverse iteration on MATRIX A (see README)!\n\n");
  double eigVal = 54; gsl_vector_memcpy(x,tmpx);
  steps = inv_iter(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 2, eigvec = [1 1 1]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);

  eigVal = 0.1; gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
  steps = inv_iter(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 0.5, eigvec = [1 0 2]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);
  
  eigVal = -10; gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
  steps = inv_iter(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = -0.5, eigvec = [-1 2 0]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);
  
  end_t = clock();
  double diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> timing: %f ms <<<<\n\n", diff_t*1000);

  start_t = clock();

  printf("Testing inverse iteration on MATRIX B (see README)!\n\n");
  eigVal = 120; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = inv_iter(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 4, eigvec = [9 4 6]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);
  
  eigVal = 2; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = inv_iter(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 3, eigvec = [1 0 1]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);
  
  eigVal = -120; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = inv_iter(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 0, eigvec = [-1 0 2]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);

  end_t = clock();
  diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> timing: %f ms <<<<\n\n", diff_t*1000);


  start_t = clock();

  
  printf("Testing inverse iteration with dynamic Rayleigh quotient on MATRIX A (see README)!\n\n");
  eigVal = 54; gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
  steps = rayleigh(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 2, eigvec = [1 1 1]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);


  eigVal = 0.1; gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
  steps = rayleigh(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 0.5, eigvec = [1 0 2]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);

  eigVal = -10; gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
  steps = rayleigh(A,x,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = -0.5, eigvec = [-1 2 0]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(x); printf("Steps = %i\n\n\n",steps);

  end_t = clock();
  diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> timing: %f ms <<<<\n\n", diff_t*1000);


  printf("Testing inverse iteration with dynamic Rayleigh quotient on MATRIX B (see README)!\n\n");
  eigVal = 120; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = rayleigh(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 4, eigvec = [9 4 6]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);

  eigVal = 2; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = rayleigh(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 3, eigvec = [1 0 1]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);

  eigVal = -120; gsl_matrix_memcpy(B, tmpB); gsl_vector_memcpy(y,tmpy);
  steps = rayleigh(B,y,eigVal,epsilon,maxIter,&eigenValue);
  printf("Expected eigval = 0, eigvec = [-1 0 2]\n");
  printf("Result eigval = %lg (guess=%lg), eigvec =\n",eigenValue, eigVal);
  printv(y); printf("Steps = %i\n\n\n",steps);

  end_t = clock();
  diff_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
  printf("\n>>>> timing: %f ms <<<<\n\n", diff_t*1000);

  /* See how number of steps increases for how bad initial guess is. */
  eigVal = 2;
  int stepiter=0;
  int steprayl=0;
  double eig1, eig2;
  fprintf(stderr, "#shift stepIter stepRayleigh eig(iter) eig(rayleigh)\n");
  for(int k=0; k<1000; k++) {
    gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
    stepiter = inv_iter(A,x,eigVal + (double)k,epsilon,maxIter,&eig1);
    gsl_matrix_memcpy(A, tmpA); gsl_vector_memcpy(x,tmpx);
    steprayl = rayleigh(A,x,eigVal + (double)k,epsilon,maxIter,&eig2);
    fprintf(stderr, "%i %i %i %lg %lg\n", k, stepiter, steprayl, eig1, eig2);
  }



  gsl_matrix_free(A);
  gsl_vector_free(x);
  gsl_matrix_free(tmpA);
  gsl_vector_free(tmpx);
  gsl_matrix_free(B);
  gsl_vector_free(y);
  gsl_matrix_free(tmpB);
  gsl_vector_free(tmpy);
  return 0;
}
