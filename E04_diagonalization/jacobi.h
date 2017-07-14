
void restore_matrix     (gsl_matrix* A);
int jacobi_sweep        (gsl_matrix* A, gsl_vector* eig, gsl_matrix* V);
int jacobi_sweep_row    (gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh);
int jacobi_sweep_row_max(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int nEig, int LowOrHigh);

void printm(const gsl_matrix* A);
void printv(const gsl_vector* v);
void matrix_set_rand(gsl_matrix* A);
void vector_set_rand(gsl_vector *v);
