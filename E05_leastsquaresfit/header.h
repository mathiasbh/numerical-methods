void   qr_givens_decomp(gsl_matrix* A);
void   qr_givens_solve (gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
double qr_givens_absdet(const gsl_matrix* QR);
void   qr_givens_invert(const gsl_matrix* QR, gsl_matrix* Ainverse); 
void   qr_givens_unpackQ(gsl_matrix* QR, gsl_matrix *Q);
void   qr_givens_unpackR(gsl_matrix* QR, gsl_matrix *R);

void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, 
	   int nf, double f(int k, double x), gsl_vector* c, gsl_matrix* S);

double func(int i, double x);


void restore_matrix     (gsl_matrix* A);
int jacobi_sweep        (gsl_matrix* A, gsl_vector* eig, gsl_matrix* V);
int jacobi_sweep_row    (gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh);
int jacobi_sweep_row_max(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh);


void lssvd(gsl_vector* x, gsl_vector* y, gsl_vector* dy, 
	   int nf, double f(int k, double x), gsl_vector* c, gsl_matrix* Sigma);
