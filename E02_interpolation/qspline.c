/* Quadratic spline
   S_i(x) = yi + bi(x-xi) + ci(x-xi)**2
   bi = pi - ci*dxi */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int binary_search(int n, double *x, double z);
typedef struct {int n; double *x, *y, *b, *c;} qspline;


qspline* qspline_alloc(int n, double* x, double* y) {
  /* Allocate memory to structure */
  qspline *s = (qspline*)malloc(sizeof(qspline));
  s->x = (double*)malloc(n*sizeof(double)); 
  s->y = (double*)malloc(n*sizeof(double));
  s->b = (double*)malloc((n-1)*sizeof(double));
  s->c = (double*)malloc((n-1)*sizeof(double));
  s->n = n;

  /* Set qspline structure "x" and "y" to the "x" and "y" from argument */
  for(int i=0; i<n; i++) {
    s->x[i] = x[i];
    s->y[i] = y[i];
  }

  /* Calculate p = dy/dx */
  double h[n-1],  p[n-1];
  for(int i=0; i<n-1; i++) {
    h[i] = x[i+1] - x[i];
    p[i] = (y[i+1] - y[i])/h[i];
  }
  
  /* Calculate c with recursion up */
  s->c[0] = 0;
  for(int i=0; i<n-2; i++) { s->c[i+1] = (p[i+1] - p[i] - s->c[i]*h[i])/h[i+1]; }

  s->c[n-2]/= 2; /* last index divided by two to get average */
  /* Calculate c with recursion down */
  for(int i=n-3; i>=0; i--) { s->c[i] = (p[i+1] - p[i] - s->c[i+1]*h[i+1])/h[i]; }

  /* Calculate b[i] = p[i] -  c[i]*h[i] */
  for(int i=0; i<n-1; i++) { s->b[i] = p[i] - s->c[i]*h[i]; }
  return s;
}

double qspline_evaluate(qspline *s, double z) {
  /* Si(x) = yi+bi(x-xi)+ci(x-xi)^2 */
  assert(z >= s->x[0] && z <= s->x[s->n-1]);

  /* Binary search */
  int i = binary_search(s->n, s->x, z);

  double h = z-s->x[i];
  return s->y[i] + s->b[i]*h + s->c[i]*h*h;
}

double qspline_derivative(qspline *s, double z) {
  /* Evaluate derivative at x=z */
  assert(z >= s->x[0] && z <= s->x[s->n-1]);

  /* Binary search */
  int i = binary_search(s->n, s->x, z);
  
  double h = z - s->x[i];
  return s->b[i] + 2*h*s->c[i];
}

double qspline_integral(qspline *s, double z) {
  /* Evaluate integral from x[0] to x=z */
  assert(z>=s->x[0] && z<=s->x[s->n-1]);
  double IntegralSum = 0, h; 

  /* Binary search */
  int i = binary_search(s->n, s->x, z);

  for(int j=0; j<i; j++) {
    h = s->x[j+1] - s->x[j];
    IntegralSum += s->y[j]*h + s->b[j]*h*h/2.0 + s->c[j]*h*h*h/3.0;
  }
  
  h = z - s->x[i];
  IntegralSum += s->y[i]*h + s->b[i]*h*h/2.0 + s->c[i]*h*h*h/3.0;
  return IntegralSum;
}

void qspline_free(qspline *s) {
  free(s->x);
  free(s->y);
  free(s->b);
  free(s->c);
  free(s);
}
