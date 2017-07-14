
void newton_minimize(void f(gsl_vector* x, double* fx),
		     void grad(gsl_vector* x, gsl_vector* df),
		     void hess(gsl_vector* x, gsl_matrix* H),
		     gsl_vector* x, double dx, double epsilon);

void newton_minimize_quasi(void f(gsl_vector* x, double* fx),
			   void grad(gsl_vector* x, gsl_vector* df),
			   gsl_vector* x, double dx, double epsilon);

void fRosenbrock(gsl_vector* w, double* fw);
void fRosenbrock_grad(gsl_vector* w, gsl_vector* fw);
void fRosenbrock_hess(gsl_vector* w, gsl_matrix* H);

void fHimmelblau(gsl_vector* w, double* fw);
void fHimmelblau_grad(gsl_vector* w, gsl_vector* fw);
void fHimmelblau_hess(gsl_vector* w, gsl_matrix* H);

gsl_matrix* jHimmelBlau_grad(gsl_vector* w);
gsl_matrix* jRosenbrock_grad(gsl_vector* w);
void fBeale(gsl_vector* w, double* fw);
void fBeale_grad(gsl_vector* w, gsl_vector* fw);
