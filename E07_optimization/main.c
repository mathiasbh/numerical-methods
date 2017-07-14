
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "header.h"

int nFunctionCalls;

int main() {
  gsl_vector* x = gsl_vector_alloc(2);
  double dx = 1e-7, epsilon = 1e-8;
  /* Exercise A: minimum of Rosenbrock (x,y) */
  /* Guess */
  gsl_vector_set(x, 0, 0.0);
  gsl_vector_set(x, 1, 0.0);
  
  printf("--- Minimization of Rosenbrock's valley ---:\n");
  printf("Part A: user specified 1st and 2nd derivative\n");
  nFunctionCalls = 0;
  newton_minimize(fRosenbrock, fRosenbrock_grad, fRosenbrock_hess, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions (x,y) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise B: minimum of Rosenbrock (x,y) */
  gsl_vector_set(x, 0, 0.0);
  gsl_vector_set(x, 1, 0.0);
  
  printf("Part B: inverse hessian\n");
  nFunctionCalls = 0;
  newton_minimize_quasi(fRosenbrock, fRosenbrock_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions (x,y) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));


  /* Exercise A: Minimum of Himmelblau (x,y) */
  gsl_vector_set(x, 0, 4.0);
  gsl_vector_set(x, 1, 4.0);
  
  printf("--- Minimization of Himmelblau's function ---:\n");
  printf("Part A: user specified 1st and 2nd derivative\n");
  nFunctionCalls = 0;
  newton_minimize(fHimmelblau, fHimmelblau_grad, fHimmelblau_hess, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions (x,y) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise B: minimum of Himmelblau (x,y) */
  gsl_vector_set(x, 0, 4.0);
  gsl_vector_set(x, 1, 4.0);
  
  printf("Part B: inverse hessian\n");
  nFunctionCalls = 0;
  newton_minimize_quasi(fHimmelblau, fHimmelblau_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions (x,y) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise B: Find minimum of Beale's function // check effectiveness if root finding */
  gsl_vector_set(x, 0, 4.0); gsl_vector_set(x, 1, 0.0);
  printf("--- Minimum of Beale's function (other interesting example) ---:\n");
  printf("Part B\n"); nFunctionCalls = 0;
  newton_minimize_quasi(fBeale, fBeale_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  printf("Comparing root finding method see ROOT_FIND_COMPARISON.txt.\nEquivalent initial guesses are used.\n");
  
  /* double t[] = {0.23,1.29,2.35,3.41,4.47,5.53,6.59,7.65,8.71,9.77}; */
  /* double y[] = {4.64,3.38,3.01,2.55,2.29,1.67,1.59,1.69,1.38,1.46}; */
  /* double e[] = {0.42,0.37,0.34,0.31,0.29,0.27,0.26,0.25,0.24,0.24}; */
  /* int N = sizeof(t)/sizeof(t[0]); */
  
  gsl_vector_free(x);
  return 0;
}
