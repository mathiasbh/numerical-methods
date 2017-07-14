#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void   qr_givens_decomp(gsl_matrix* A);
void   qr_givens_solve (gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
double qr_givens_absdet(const gsl_matrix* QR);
void   qr_givens_invert(const gsl_matrix* QR, gsl_matrix* Ainverse); 
void   qr_givens_unpackQ(gsl_matrix* QR, gsl_matrix *Q);
void   qr_givens_unpackR(gsl_matrix* QR, gsl_matrix *R);


void printm(const gsl_matrix* A);
void printv(const gsl_vector* v);
void matrix_set_rand(gsl_matrix* A);
void vector_set_rand(gsl_vector *v);


int inv_iter(gsl_matrix* A, gsl_vector* b, double s, double eps, int maxIter, double* eigenValue);
int rayleigh(gsl_matrix* A, gsl_vector* b, double s, double eps, int maxIter, double* eigenValue);
