
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>


void qr_givens_decomp(gsl_matrix* A) { /* QR -> A */
  int n = A->size1;
  int m = A->size2;
  for(int q=0; q<m; q++) { /* Run over all off diagonal in lower triangle of A, column per column */
    for(int p=q+1; p<n; p++) {
      /* tan theta = xp/xq <=> theta = atan(xp/xq) */
      double theta = atan2(gsl_matrix_get(A, p, q), gsl_matrix_get(A, q, q));
      for(int k=q; k<m; k++) {
	double xq = gsl_matrix_get(A, q, k);
	double xp = gsl_matrix_get(A, p, k);
	
	gsl_matrix_set(A, q, k, xq*cos(theta) + xp*sin(theta));
	gsl_matrix_set(A, p, k,-xq*sin(theta) + xp*cos(theta));
      }
      /* Store angle in empty lower triangle */
      gsl_matrix_set(A, p, q, theta);
    }
  }
}

void qr_givens_QTvec(gsl_matrix* QR, gsl_vector* v) {
 
  int n = QR->size1;
  int m = QR->size2;

  /* Calculate Gb by successively applying Given rotations - store in v */
  for(int q=0; q<m; q++) {
    for(int p=q+1; p<n; p++) {
      double theta = gsl_matrix_get(QR, p, q);
      double vq = gsl_vector_get(v, q);
      double vp = gsl_vector_get(v, p);
      gsl_vector_set(v, q, vq*cos(theta) + vp*sin(theta));
      gsl_vector_set(v, p,-vq*sin(theta) + vp*cos(theta));
    }
  }
}
 
void qr_givens_solve(gsl_matrix* QR, gsl_vector* b, gsl_vector* x) {
  /* Solves Ax=b - stores result in vector x */
  /* Use Givens angles (stored in QR) to make equation Rx = Gb */ 
  /* Solve Rx = Gb by backsubstitution - store result in x */

  int m = QR->size2;
  qr_givens_QTvec(QR, b);

  for(int i=m-1; i>=0; i--) {
    double SumUy = 0;
    for(int k=i+1; k<m; k++) {
      SumUy -= gsl_matrix_get(QR, i, k) * gsl_vector_get(x, k);
    }
    gsl_vector_set(x, i, (gsl_vector_get(b,i) + SumUy)/gsl_matrix_get(QR, i, i));
  }
}

void gsl_vector_set_unit(gsl_vector* e, int k) {
  for(int i=0; i<e->size; i++) gsl_vector_set(e,i,0);
  gsl_vector_set(e,k,1);
}

void qr_givens_unpackQ(gsl_matrix* QR, gsl_matrix* Q) {
  gsl_vector* ei = gsl_vector_alloc(QR->size1);
  for(int i=0; i<QR->size1; i++) {
    gsl_vector_set_unit(ei, i);
    qr_givens_QTvec(QR, ei);
    for(int j=0; j<QR->size2; j++) {
      gsl_matrix_set(Q, i, j, gsl_vector_get(ei, j));
    }
  }
  gsl_vector_free(ei);
}

void qr_givens_unpackR(gsl_matrix* QR, gsl_matrix* R) {
  for(int i=0; i<R->size2; i++) {
    for(int j=0; j<R->size2; j++) {
      if(j<=i) gsl_matrix_set(R,j,i, gsl_matrix_get(QR,j,i));
      else gsl_matrix_set(R,j,i,0);
    }
  }
}

double qr_givens_absdet(const gsl_matrix* QR) {
  /* Calculate determinant */
  int n = QR->size1;
  double ProductSum = 1;
  for(int i=0; i<n; i++) {
    ProductSum *= gsl_matrix_get(QR, i, i);
  }
  
  return fabs(ProductSum);
}

void qr_givens_invert(gsl_matrix* QR, gsl_matrix* Ainverse) {
  int n = Ainverse->size1;

  gsl_matrix_set_identity(Ainverse); /* Set Ainverse to identity matrix */
  gsl_vector* x = gsl_vector_calloc(n);

  for(int i=0; i<n; i++) {
    gsl_vector_view ei = gsl_matrix_column(Ainverse,i);
    qr_givens_solve(QR, &ei.vector, x);
    gsl_matrix_set_col(Ainverse, i, x);    
  }

  gsl_vector_free(x);
}
