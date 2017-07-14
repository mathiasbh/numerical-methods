#define RND ((double)rand()/(double)RAND_MAX)
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void printm(const gsl_matrix* A) {
  int n = A->size1; /* Row */
  int m = A->size2; /* Column */

  for(int i=0; i<n; i++) {
    for(int j=0; j<m; j++) {
      printf("%10.3f ", gsl_matrix_get(A,i,j)); /* Print all j's in row i */
    } printf("\n");
  } printf("\n");
}

void printv(const gsl_vector* v) {
  int n = v->size;

  for(int i=0; i<n; i++) {
    printf("%10.3f\n", gsl_vector_get(v, i));
  } printf("\n");
}

void matrix_set_rand(gsl_matrix* A) {
  for(int i=0; i<A->size1; i++) {
    gsl_matrix_set(A, i, i, RND);
    for(int j=i+1; j<A->size2; j++) {
      gsl_matrix_set(A, i, j, RND);
      gsl_matrix_set(A, j, i, gsl_matrix_get(A, i, j));
    }
  }
}

void vector_set_rand(gsl_vector *v) {
  for(int i=0; i<v->size; i++) {
    gsl_vector_set(v, i, RND);
  }
}
