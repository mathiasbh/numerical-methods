#include <math.h>
#include <assert.h>
#include <stdio.h>
#define M_PI 3.14159265358979323846264338

double adaptOpen(double f(double), double a, double b,
	       double acc, double eps, double f2, double f3,
	       double nrec, double* err) {
  assert(nrec < 1e7);
  double f1 = f(a+(b-a)/6.0);   // equ 48, f2,f3 need not be recalculated
  double f4 = f(a+5*(b-a)/6.0); // they are inherited from previous step
  double Q  = (2*f1+f2+f3+2*f4)/6.0 * (b-a); // eqn 49 and 51
  double q  = (f1+f2+f3+f4)/4.0 * (b-a);
  double tol= acc + eps*fabs(Q); // eqn 47
  *err = fabs(Q-q); // eqn 46
  if (*err < tol) {
    return Q;
  } else {
    double err1, err2;
    double Q1=adaptOpen(f,a,(a+b)/2.0,acc/sqrt(2.0),eps,f1,f2,nrec+1,&err1);
    double Q2=adaptOpen(f,(a+b)/2.0,b,acc/sqrt(2.0),eps,f3,f4,nrec+1,&err2);
    *err = sqrt(err1*err1 + err2*err2);
    return Q1+Q2;
  }
}

double intOpen(double f(double), double a, double b,
	     double acc, double eps, double* err) {
  double f2,f3;
  int nrec = 0;
  if (isinf(a) && isinf(b)) { // -inf to +inf
    double ft(double t) {return (f((1.0-t)/t) + f((t-1.0)/t))/t/t;} /* Eqn 58 */
    a = 0; b = 1;
    f2 = ft(a+2*(b-a)/6.0);
    f3 = ft(a+4*(b-a)/6.0);
    return adaptOpen(ft,a,b,acc,eps,f2,f3,nrec,err);
  } else if (isinf(b)) { // a to inf
    double ft(double t) {return f(a+(1.0-t)/t)/(t*t);}  /* Eqn 60 */
    a=0;
    b=1;
    f2 = ft(a+2*(b-a)/6.0);
    f3 = ft(a+4*(b-a)/6.0);
    return adaptOpen(ft,0,b,acc,eps,f2,f3,nrec,err);
  } else if (isinf(a)) { // -inf to a
    double ft(double t) {return f(a-(1.0-t)/t)/(t*t);} /* Eqn 62 */
    a = 0; b = 1;
    f2 = ft(a+2*(b-a)/6.0);
    f3 = ft(a+4*(b-a)/6.0);
    return adaptOpen(ft,a,b,acc,eps,f2,f3,nrec,err);
  } else {
    f2 = f(a+2*(b-a)/6.0); // eqn 48
    f3 = f(a+4*(b-a)/6.0);
    return adaptOpen(f,a,b,acc,eps,f2,f3,nrec,err);
  }
}

double intOpen2d(int size, double f(double* x), double* a, double* b, double acc, double eps, double* error) {
  double a0 = a[0], a1 = a[1];
  double b0 = b[0], b1 = b[1];
  double F(double x) {
    double fph(double y) {
      double vec[] = {x, y};
      return f(vec);
    }
    return intOpen(fph,a1,b1,acc,eps,error);
  }
  return intOpen(F,a0,b0,acc,eps,error);
}
