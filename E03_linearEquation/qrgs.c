/* QR-decomposition by modified Gram-Schmidt orthogonalization
   A = QR
   Convert linear system Ax=b into trinagular form Rx = Q^T*b
   Q: orthogonal matrix, Q^T*Q = 1
   R: triangular matrix */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {
  /* Perform in-place modified Gram-Schmidt orthogonalization of 
     matrix A (turn A into Q) and fill matrix R */

  int m = A->size2;
  for(int i=0; i<m; i++) { 
    gsl_vector_view ai = gsl_matrix_column(A, i); /* Create vector view */
    /* Calculate norm ||x||_2 = sqrt(sum(x_i^2)) */
    double NormVector = gsl_blas_dnrm2(&ai.vector);
    
    /* R_ii = NormVector */
    gsl_matrix_set(R, i, i, NormVector);

    /* q_i = a_i / R_ii */
    gsl_vector_scale(&ai.vector, 1/NormVector);

    for(int j=i+1; j<m; j++) { /* Off diagonal */
      gsl_vector_view aj = gsl_matrix_column(A, j);

      /* R_ij = q_i^T * a_j */
      double ScalerProduct = 0;
      gsl_blas_ddot(&ai.vector, &aj.vector, &ScalerProduct);

      /* a_j = a_j - q_i * R_ij */
      gsl_blas_daxpy( -ScalerProduct, &ai.vector, &aj.vector); /* aj = -scalarP*ai + aj */

      gsl_matrix_set(R, i, j, ScalerProduct);
    }
  }
}


void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, const gsl_vector* b, gsl_vector* x) {
  /* Given Q and R, solve QRx=b, by applying Q^T on righthand side 
     of b and perform back-substitution */
  /* Rx = Q^T*b */
  int n = R->size2;
  gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x); /* Q^T*b, sum y = alpha*Q^T*b + beta*y */
  for(int i=n-1; i>=0; i--) { /* Back substitution */
    double SumUy = 0;
    for(int k=i+1; k<n; k++) {
      SumUy -= gsl_matrix_get(R, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, (gsl_vector_get(x,i) + SumUy)/gsl_matrix_get(R, i, i));
  }
  
}

double qr_gs_absdet(const gsl_matrix *R) {
  /* Given R, calculate absolute value of determinant of matrix A */
  /* det A = det Q det R
     det Q = 1, since it is orthogonal 
     det A = det R = productSum R_ii */

  int m = R->size2;
  double ProductSum = 1;
  for(int i=0; i<m; i++) {
    ProductSum *= gsl_matrix_get(R, i, i);
  }
  
  return fabs(ProductSum);
}

void qr_gs_invert(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* Ainverse) {
  /* Given Q,R calculate inverse of matrix A into Ainverse */
  /* Inverse A^-1 can be calculated by solving n linear equations
     A x_i = unitvector_i | i=1,...,n */
    
  int n = Ainverse->size1;

  gsl_vector* b = gsl_vector_calloc(n);
  gsl_vector* x = gsl_vector_calloc(n);

  /* Set Ainverse to identity matrix */
  for(int i=0; i<n; i++) {
    gsl_vector_set(b, i, 1.0); /* Make b into a unit vector */
    qr_gs_solve(Q, R, b, x); /* Solve equation */
    gsl_vector_set(b, i, 0.0); /* Set b to zero matrix again for next iteration */
    gsl_matrix_set_col(Ainverse, i, x);    
  }
  gsl_vector_free(b);
  gsl_vector_free(x);
}
