/* Cubic spline
   S_i(x) = yi + bi(x-xi) + ci(x-xi)**2 + di(x-xi)**3
   bi = pi - ci*dxi */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int binary_search(int n, double *x, double z);
typedef struct {int n; double *x, *y, *b, *c, *d;} cspline;


cspline* cspline_alloc(int n, double* x, double* y) {
  /* Allocate memory to structure */
  cspline *s = (cspline*)malloc(sizeof(cspline));
  s->x = (double*)malloc(n*sizeof(double)); 
  s->y = (double*)malloc(n*sizeof(double));
  s->b = (double*)malloc(n*sizeof(double));
  s->c = (double*)malloc((n-1)*sizeof(double));
  s->d = (double*)malloc((n-1)*sizeof(double));
  s->n = n;

  /* Set cspline structure "x" and "y" to the "x" and "y" from argument */
  for(int i=0; i<n; i++) {
    s->x[i] = x[i];
    s->y[i] = y[i];
  }

  /* Calculate p = dy/dx */
  double h[n-1],  p[n-1];
  for(int i=0; i<n-1; i++) {
    h[i] = x[i+1] - x[i]; assert(h[i] > 0);
    p[i] = (y[i+1] - y[i])/h[i];
  }
  
  /* Build tridiagonal system */
  /* | D1 Q1 0  0    ... | b1 |   | B1 | */
  /* | 1  D2 Q2 0  0 ... | .. |   | .. | */
  /* | 0  1  D3 Q3 0 ... | .. | = | .. | */
  /* | ...               | .. |   | .. | */
  /* | ... ...  0  1 Dn  | bn |   | Bn | */
  double D[n], Q[n-1], B[n];
  for(int i=0; i<n-2; i++) {
    D[i+1] = 2*h[i]/h[i+1]+2;
    Q[i+1] = h[i]/h[i+1];
    B[i+1] = 3*(p[i] + p[i+1]*h[i]/h[i+1]);
  }
  D[0] = 2;
  Q[0] = 1;
  B[0] = 3*p[0];
  D[n-1] = 2;
  B[n-1] = 3*p[n-2];
  
  /* Perform Gaussian elimination */
  for(int i=1; i<n; i++) {
    D[i] -= Q[i-1]/D[i-1]; /* D1_tilde = D1, Di_tilde = Di - Q(i-1)/D(i-1)_tilde */
    B[i] -= B[i-1]/D[i-1]; /* B1_tilde = B1, Bi_tilde = Bi - B(i-1)_tilde/D(i-1)_tilde */
  }
 
  /* Calculate b, c and d using back substitution */
  s->b[n-1] = B[n-1]/D[n-1];
  for(int i=n-2; i>=0; i--) s->b[i] = (B[i]-Q[i]*s->b[i+1])/D[i];
  for(int i=0; i<n-1; i++) {
    s->c[i] = (-2*s->b[i] - s->b[i+1] + 3*p[i])/h[i];
    s->d[i] = (s->b[i] + s->b[i+1] - 2*p[i])/h[i]/h[i];
  }
  return s;
}

double cspline_evaluate(cspline *s, double z) {
  /* Si(x) = yi+bi(x-xi)+ci(x-xi)^2+di(x-xi)^3 */
  assert(z >= s->x[0] && z <= s->x[s->n-1]);

  /* Binary search */
  int i = binary_search(s->n, s->x, z);

  double h = z-s->x[i];
  return s->y[i] + s->b[i]*h + s->c[i]*h*h + s->d[i]*h*h*h;
}

double cspline_derivative(cspline *s, double z) {
  /* Evaluate derivative at x=z */
  assert(z >= s->x[0] && z <= s->x[s->n-1]);

  /* Binary search */
  int i = binary_search(s->n, s->x, z);
  
  double h = z - s->x[i];
  return s->b[i] + 2*h*s->c[i] + 3*h*h*s->d[i];
}

double cspline_integral(cspline *s, double z) {
  /* Evaluate integral from x[0] to x=z */
  assert(z >= s->x[0] && z <= s->x[s->n-1]);

  double IntegralSum = 0, h;

  /* Binary search */
  int i = binary_search(s->n, s->x, z);
  
  for(int j=0; j<i; j++) {
    h = s->x[j+1] - s->x[j];
    IntegralSum += s->y[j]*h + s->b[j]*h*h/2 + s->c[j]*h*h*h/3 + s->d[j]*h*h*h*h/4;
  }
  
  /* Integral at x=z */
  h = z - s->x[i];
  IntegralSum += s->y[i]*h + s->b[i]*h*h/2 + s->c[i]*h*h*h/3 + s->d[i]*h*h*h*h/4;

  return IntegralSum;
}

void cspline_free(cspline *s) {
  free(s->x);
  free(s->y);
  free(s->b);
  free(s->c);
  free(s->d);
  free(s);
}
