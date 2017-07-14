/* Matrix diagonalization using Jacobi */
/* For real symmetric matrices */

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#define M_PI 3.14159265358979323846264338


void restore_matrix(gsl_matrix* A) {
  int n = A->size1;
  for(int i=0; i<n; i++) {
    for(int j=i+1; j<n; j++) {
      gsl_matrix_set(A, i, j, gsl_matrix_get(A, j, i));
    }
  }
}

int jacobi_sweep(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V) {
  /* Sweep: sequence of Jacobi rotations applted to all non-diagonal elements. */
  /* Eigenvectors calculated as V = 1*J1*J2..., J: Jacobi matrices */
  /* A -> A' = J^T*A*J, A'pq = A'qp = 0 */
  /* Upper triangle of A is destroyed */
  /* eig and V accumulate eigenvalues and vectors */

  assert(A->size1 == A->size2);

  int changed; /* Check if elements have changed */
  int rotations = 0; /* Number of rotations */
  int n = A->size1;
  int p, q;
  
  /* V is initially the identity */
  gsl_matrix_set_identity(V); 

  /* Diagonal of A is copied into eig */
  for(int i=0; i<n; i++) {
    gsl_vector_set(eig, i, gsl_matrix_get(A, i, i));
  }

  do {
    changed = 0;
    for(p=0; p<n; p++) {
      for(q=p+1; q<n; q++) { /* Upper off diagonal: p is row q is column */
	double App=gsl_vector_get(eig,p);
	double Aqq=gsl_vector_get(eig,q);
	double Apq=gsl_matrix_get(A,p,q);
	double phi=0.5*atan2(2*Apq,Aqq-App); /* => Apq' = 0 */
	double c = cos(phi);
	double s = sin(phi);

	double Appdash = c*c*App - 2*s*c*Apq + s*s*Aqq;
	double Aqqdash = s*s*App + 2*s*c*Apq + c*c*Aqq;

	/* Check if elements has changed. If changed, apply Jacobi 
	   since eigenvalues have yet to converge */
	if(Appdash != App || Aqqdash != Aqq) {
	  changed = 1; 
	  rotations++;
	  gsl_vector_set(eig, p, Appdash);
	  gsl_vector_set(eig, q, Aqqdash);
	  gsl_matrix_set(A, p, q, 0.0); /* p,q element equals zero */
	


	  for(int i=0; i<p; i++){
	    double Aip = gsl_matrix_get(A, i, p);
	    double Aiq = gsl_matrix_get(A, i, q);
	    gsl_matrix_set(A, i, p, c*Aip - s*Aiq);
	    gsl_matrix_set(A, i, q, c*Aiq + s*Aip);
	  }

	  for(int i=p+1; i<q; i++){
	    double Api = gsl_matrix_get(A, p, i);
	    double Aiq = gsl_matrix_get(A, i, q);
	    gsl_matrix_set(A, p, i, c*Api - s*Aiq);
	    gsl_matrix_set(A, i, q, c*Aiq + s*Api);
	  }

	  for(int i=q+1; i<n; i++){ /*  */
	    double Api = gsl_matrix_get(A, p, i);
	    double Aqi = gsl_matrix_get(A, q, i);
	    gsl_matrix_set(A, p, i, c*Api - s*Aqi);
	    gsl_matrix_set(A, q, i, c*Aqi + s*Api);
	  }

	  for(int i=0; i<n; i++){
	    double vip = gsl_matrix_get(V, i, p);
	    double viq = gsl_matrix_get(V, i, q);
	    gsl_matrix_set(V, i, p, c*vip - s*viq);
	    gsl_matrix_set(V, i, q, c*viq + s*vip);
	  }
	}
      }
    }
  } while (changed != 0);
  return rotations;
}


/* Equivalent as jacobi_sweep(), but by running over a row at a time */
/* nEig is the number of eigenvalues to be found */
/* LowOrHigh >= 0 => find largest eigenvalues
   LowOrHigh <  0 => find lowest eigenvalues */

int jacobi_sweep_row(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh) {
  assert(A->size1 == A->size2);
  assert(A->size1 > nEig && nEig > 0);

  int changed; /* Check if elements have changed */
  int rotations = 0; /* Number of rotations */
  int n = A->size1;
  int p, q;
  
  /* V is initially the identity */
  gsl_matrix_set_identity(V); 

  /* Diagonal of A is copied into eig */
  for(int i=0; i<n; i++) {
    gsl_vector_set(eig, i, gsl_matrix_get(A, i, i));
  }

  for(p=0; p<nEig; p++) { /* Run over each row */
    do {
      changed = 0;
      
      for(q=p+1; q<n; q++) { /* Upper off diagonal: p is row q is column */
	double App=gsl_vector_get(eig,p);
	double Aqq=gsl_vector_get(eig,q);
	double Apq=gsl_matrix_get(A,p,q);
	double phi; /* => Apq' = 0 */
	if(LowOrHigh < 0) phi=0.5*atan2(2*Apq,Aqq-App); /* Low eigenvalues */
	else phi=0.5*atan2(-2*Apq,App-Aqq);             /* High eigenvalues */
	double c = cos(phi);
	double s = sin(phi);

	double Appdash = c*c*App - 2*s*c*Apq + s*s*Aqq;
	double Aqqdash = s*s*App + 2*s*c*Apq + c*c*Aqq;

	/* Check if elements has changed. If changed, apply Jacobi 
	   since eigenvalues have yet to converge */
	if(Appdash != App || Aqqdash != Aqq) {
	  changed = 1; 
	  rotations++;
	  gsl_vector_set(eig, p, Appdash);
	  gsl_vector_set(eig, q, Aqqdash);
	  gsl_matrix_set(A, p, q, 0.0); /* p,q element equals zero */

	  for(int i=p+1; i<q; i++){
	    double Api = gsl_matrix_get(A, p, i);
	    double Aiq = gsl_matrix_get(A, i, q);
	    gsl_matrix_set(A, p, i, c*Api - s*Aiq);
	    gsl_matrix_set(A, i, q, c*Aiq + s*Api);
	  }

	  for(int i=q+1; i<n; i++){ /*  */
	    double Api = gsl_matrix_get(A, p, i);
	    double Aqi = gsl_matrix_get(A, q, i);
	    gsl_matrix_set(A, p, i, c*Api - s*Aqi);
	    gsl_matrix_set(A, q, i, c*Aqi + s*Api);
	  }

	  for(int i=0; i<n; i++){
	    double vip = gsl_matrix_get(V, i, p);
	    double viq = gsl_matrix_get(V, i, q);
	    gsl_matrix_set(V, i, p, c*vip - s*viq);
	    gsl_matrix_set(V, i, q, c*viq + s*vip);
	  }
	}
      }
    } while (changed != 0);
  }
  return rotations;
}



/* Equivalent as jacobi_sweep(), but by running over a row at a time 
   finding the largest first and eliminating it.
   nEig is the number of eigenvalues to be found */

int jacobi_sweep_row_max(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh) {
  assert(A->size1 == A->size2);
  assert(A->size1 > nEig && nEig > 0);

  int changed; /* Check if elements have changed */
  int rotations = 0; /* Number of rotations */
  int n = A->size1;
  int p, q;
  
  /* V is initially the identity */
  gsl_matrix_set_identity(V); 

  /* Diagonal of A is copied into eig */
  for(int i=0; i<n; i++) {
    gsl_vector_set(eig, i, gsl_matrix_get(A, i, i));
  }

  for(p=0; p<nEig; p++) { /* Run over row */
    do {
      changed = 0;
      
      /* Find largest element in row instead of running over entire row */
      q = p + 1;
      double CurrentElement = fabs(gsl_matrix_get(A,p,q)); /* Test if element is larger than others */
      for(int i=p+2; i<n; i++) {
	double TestElement = fabs(gsl_matrix_get(A,p,i));
	if(CurrentElement < TestElement) {
	  q = i; /* Move to i, since this is current largest number */
	  CurrentElement = TestElement;
	}
      }
	
      
      double App=gsl_vector_get(eig,p);
      double Aqq=gsl_vector_get(eig,q);
      double Apq=gsl_matrix_get(A,p,q);

      double phi=0.5*atan2(2*Apq,Aqq-App);

      if(LowOrHigh > 0) phi += M_PI/2;     /* high eigenvalues, low is automatic */
      

      double c = cos(phi);
      double s = sin(phi);

      double Appdash = c*c*App - 2*s*c*Apq + s*s*Aqq;
      double Aqqdash = s*s*App + 2*s*c*Apq + c*c*Aqq;

      /* Check if elements has changed. If changed, apply Jacobi 
	 since eigenvalues have yet to converge */
      if(Appdash != App || Aqqdash != Aqq) {
	changed = 1; 
	rotations++;
	gsl_vector_set(eig, p, Appdash);
	gsl_vector_set(eig, q, Aqqdash);
	gsl_matrix_set(A, p, q, 0.0); /* p,q element equals zero */

	for(int i=p+1; i<q; i++){
	  double Api = gsl_matrix_get(A, p, i);
	  double Aiq = gsl_matrix_get(A, i, q);
	  gsl_matrix_set(A, p, i, c*Api - s*Aiq);
	  gsl_matrix_set(A, i, q, c*Aiq + s*Api);
	}

	for(int i=q+1; i<n; i++){ /*  */
	  double Api = gsl_matrix_get(A, p, i);
	  double Aqi = gsl_matrix_get(A, q, i);
	  gsl_matrix_set(A, p, i, c*Api - s*Aqi);
	  gsl_matrix_set(A, q, i, c*Aqi + s*Api);
	}

	for(int i=0; i<n; i++){
	  double vip = gsl_matrix_get(V, i, p);
	  double viq = gsl_matrix_get(V, i, q);
	  gsl_matrix_set(V, i, p, c*vip - s*viq);
	  gsl_matrix_set(V, i, q, c*viq + s*vip);
	}
	
      }
    } while (changed != 0);
  }
  return rotations;
}

