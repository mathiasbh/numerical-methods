
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "newton.h"

int nFunctionCalls;

int main() {
  /* Exercise A: solve system of eqn (x,y) */
  gsl_vector* x = gsl_vector_alloc(2);
  /* Guess */
  gsl_vector_set(x, 0, 1.0);
  gsl_vector_set(x, 1, 2.0);
  
  double dx = 1e-10, epsilon = 1e-10;
  printf("--- Solve system of equations ---:\n");
  printf("Part A\n");
  nFunctionCalls = 0;
  newton_rootfind(fSystemEqn, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions (x,y) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise A: Find minimum of Beale's function */
  gsl_vector_set(x, 0, 4.0); gsl_vector_set(x, 1, 0.0);
  printf("--- Minimum of Beale's function (other interesting example) ---:\n");
  printf("Part A\n"); nFunctionCalls = 0;
  newton_rootfind(fBeale_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise A: Find minimum of Rosenbrocks valley function */
  gsl_vector_set(x, 0, 0);
  gsl_vector_set(x, 1, 0);

  /* fprintf(stderr, "Minimum of Rosenbrock's valley:\n");  */
  printf("--- Minimum of Rosenbrock's valley ---:\n");
  printf("Part A\n");
  nFunctionCalls = 0;
  newton_rootfind(fRosenbrock_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise B: Rosenbrock minima with user specified jacobian */
  gsl_vector_set(x, 0, 0);
  gsl_vector_set(x, 1, 0);

  printf("Part B\n");
  nFunctionCalls = 0;
  newton_rootfind_jacobi(fRosenbrock_grad, jRosenbrock_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise C: Rosenbrock minima with user specified jacobian */
  gsl_vector_set(x, 0, 0);
  gsl_vector_set(x, 1, 0);

  printf("Part C\n");
  nFunctionCalls = 0;
  newton_quad_rootfind(fRosenbrock_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise A: Find minimum of Himmelblau's function */
  gsl_vector_set(x, 0, 4.0); gsl_vector_set(x, 1, 4.0);
  printf("--- Minimum of Himmelblau's function ---:\n");
  printf("Part A\n"); nFunctionCalls = 0;
  newton_rootfind(fHimmelblau_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise B: Himmelblau minima with user specified jacobian */
  gsl_vector_set(x, 0, 4.0); gsl_vector_set(x, 1, 4.0);
  printf("Part B\n"); nFunctionCalls = 0;
  newton_rootfind_jacobi(fHimmelblau_grad, jHimmelBlau_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  /* Exercise C: Himmelblau minima with user specified jacobian */
  gsl_vector_set(x, 0, 4.0); gsl_vector_set(x, 1, 4.0);
  printf("Part C\n"); nFunctionCalls = 0;
  newton_quad_rootfind(fHimmelblau_grad, x, dx, epsilon);
  printf("Number of function calls: %i\n", nFunctionCalls);
  printf("Solutions minima (x1,x2) = (%lg, %lg)\n\n", gsl_vector_get(x,0), gsl_vector_get(x,1));

  gsl_vector_free(x);
  return 0;
}
