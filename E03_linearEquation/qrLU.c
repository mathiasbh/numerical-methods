/* LU-decomposition: factorization of square matrix A into product of
   lower triangular matrix L and upper U 
   A = L*U
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>


void qr_lu_decomp(gsl_matrix* A, gsl_matrix* L, gsl_matrix* U) { /* A -> LU */
  assert(gsl_matrix_get(A,0,0) != 0);

  /* LUx = b */
  int n = A->size1;
  /* int m = A->size2; */


  /* Require main diagonal of L be equal to one */
  /* Use Doolittle algoritm */
  for(int i=0; i<n; i++) {
    gsl_matrix_set(L, i, i, 1.0);

    for(int j=i; j<n; j++) {
      double Uij = gsl_matrix_get(A,i,j); /* Uij = Aij - sum Lik Ukj */

      for(int k=0; k<i; k++) {
	Uij -= gsl_matrix_get(L,i,k) * gsl_matrix_get(U,k,j);
      }
      gsl_matrix_set(U,i,j,Uij);
    }
    
    for(int j=i+1; j<n; j++) {
      double Lji = gsl_matrix_get(A,j,i); /* Lji = 1/Uii * (Aji - sum Ljk Uki) */

      for(int k=0; k<j; k++) {
	Lji -= gsl_matrix_get(L,j,k)* gsl_matrix_get(U,k,i);
      }
      /* Ensure you do not divide by 0 */
      assert(gsl_matrix_get(U,i,i) != 0);
      Lji /= gsl_matrix_get(U,i,i);
      gsl_matrix_set(L,j,i,Lji);
    }
  }

}

void qr_lu_solve(gsl_matrix* L, gsl_matrix* U, gsl_vector* b, gsl_vector* x) {
  /* Solve LUx=b by first solving Ly=b and then Ux=y */
  /* for x with two runs of forward and backward substitutions */

  int n = L->size1;

  /* First Ly=b by forward substitution */
  gsl_vector* y = gsl_vector_calloc(n);

  for(int i=0; i<n; i++) {
    double SumLy = 0;
    for(int k=0; k<i; k++) {
      SumLy -= gsl_matrix_get(L,i,k) * gsl_vector_get(y,k);
    }
    gsl_vector_set(y, i, (gsl_vector_get(b,i) + SumLy) / gsl_matrix_get(L,i,i));
  }

  /* Solve Ux=y by backsubstituion */
  for(int i=n-1; i>=0; i--) {
    double SumUx = 0;
    for(int k=i+1; k<n; k++) {
      SumUx -= gsl_matrix_get(U, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, (gsl_vector_get(y,i) + SumUx)/gsl_matrix_get(U,i,i));
  }
  
  gsl_vector_free(y);
}

void qr_lu_invert(gsl_matrix* L, gsl_matrix* U, gsl_matrix* Ainverse) {
  int n = Ainverse->size1;

  gsl_matrix_set_identity(Ainverse); /* Set Ainverse to identity matrix */
  gsl_vector* x = gsl_vector_calloc(n);

  for(int i=0; i<n; i++) {
    gsl_vector_view ei = gsl_matrix_column(Ainverse,i);
    qr_lu_solve(L, U, &ei.vector, x);
    gsl_matrix_set_col(Ainverse, i, x);    
  }

  gsl_vector_free(x);
}
