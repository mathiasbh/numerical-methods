
void qr_givens_decomp(gsl_matrix* A);
void qr_givens_QTvec(gsl_matrix* QR, gsl_vector* v);
void qr_givens_solve(gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
void gsl_vector_set_unit(gsl_vector* e, int k);
void qr_givens_unpackQ(gsl_matrix* QR, gsl_matrix* Q);
void qr_givens_unpackR(gsl_matrix* QR, gsl_matrix* R);
double qr_givens_absdet(const gsl_matrix* QR);
void qr_givens_invert(gsl_matrix* QR, gsl_matrix* Ainverse);
