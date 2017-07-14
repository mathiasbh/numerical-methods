#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
int nFunctionCalls;

void fSystemEqn(gsl_vector* w, gsl_vector* fw) {
  nFunctionCalls++;
  double A = 10000;
  double x = gsl_vector_get(w,0);
  double y = gsl_vector_get(w,1);
  
  gsl_vector_set(fw, 0, A*x*y-1.0);
  gsl_vector_set(fw, 1, exp(-x) + exp(-y) - 1.0 - 1/A);
}

void fRosenbrock_grad(gsl_vector* w, gsl_vector* fw) {
  nFunctionCalls++;
  double x = gsl_vector_get(w, 0);
  double y = gsl_vector_get(w, 1);
  
  gsl_vector_set(fw, 0, -2.0*(1.0-x) - 400.0 * (y-x*x)*x);
  gsl_vector_set(fw, 1, 200.0 * (y-x*x));
}

/* Jacobian of Rosenbrocks valley */
gsl_matrix* jRosenbrock_grad(gsl_vector* w) {
  double x = gsl_vector_get(w, 0);
  double y = gsl_vector_get(w, 1);

  gsl_matrix* jacobian = gsl_matrix_alloc(2, 2);
  gsl_matrix_set(jacobian, 0, 0, 2.0 - 400*y + 1200*x*x);
  gsl_matrix_set(jacobian, 0, 1, -400*x);
  gsl_matrix_set(jacobian, 1, 0, -400*x);
  gsl_matrix_set(jacobian, 1, 1, 200);

  return jacobian;
}


void fHimmelblau_grad(gsl_vector* w, gsl_vector* fw) {
  nFunctionCalls++;
  double x = gsl_vector_get(w, 0);
  double y = gsl_vector_get(w, 1);

  gsl_vector_set(fw, 0, 4.0*(x*x + y - 11.0)*x + 2.0*(x + y*y - 7.0));
  gsl_vector_set(fw, 1, 2.0*(x*x + y - 11.0)   + 4.0*(x + y*y - 7.0)*y);
}

/* Jacobian of Himmelblau function */
gsl_matrix* jHimmelBlau_grad(gsl_vector* w) {
  double x = gsl_vector_get(w, 0);
  double y = gsl_vector_get(w, 1);

  gsl_matrix* jacobian = gsl_matrix_alloc(2, 2);
  gsl_matrix_set(jacobian, 0, 0, 4*(3*x*x + y - 10.5));
  gsl_matrix_set(jacobian, 0, 1, 4*(x + y));
  gsl_matrix_set(jacobian, 1, 0, 4*(x + y));
  gsl_matrix_set(jacobian, 1, 1, 4*(x + 3*y*y - 6.5));

  return jacobian;
}

void fBeale_grad(gsl_vector* w, gsl_vector* fw) {
  // f(x,y) = (1.5-x+xy)^2 + (2.25-x+xy^2)^2 + (2.625-x+xy^3)^2
  // minima in f(3,0.5) = 0 according to wikipedia Test_functions_for_optimization
  nFunctionCalls++;
  double x = gsl_vector_get(w, 0);
  double y = gsl_vector_get(w, 1);
  /* Derivative of Beale's function */
  gsl_vector_set(fw, 0, 2*(y-1)*(1.5-x+x*y) + 2*(y*y-1)*(2.25-x+x*y*y) + 2*(y*y*y-1)*(2.625-x+x*y*y*y));
  gsl_vector_set(fw, 1, 2*x*(1.5-x+x*y) + 4*x*y*(2.25-x+x*y*y) + 6*x*y*y*y*(2.625-x+x*y*y*y));
}
