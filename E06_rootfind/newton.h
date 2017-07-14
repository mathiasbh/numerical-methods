
void newton_rootfind(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon);
void newton_rootfind_jacobi(void f(gsl_vector* x, gsl_vector* fx), 
			    gsl_matrix* jacobian(gsl_vector* x), gsl_vector* x, double dx, double epsilon);
void newton_quad_rootfind(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double epsilon);

void fSystemEqn(gsl_vector* w, gsl_vector* fw);
void fRosenbrock_grad(gsl_vector* w, gsl_vector* fw);
void fHimmelblau_grad(gsl_vector* w, gsl_vector* fw);
gsl_matrix* jHimmelBlau_grad(gsl_vector* w);
gsl_matrix* jRosenbrock_grad(gsl_vector* w);
void fBeale_grad(gsl_vector* w, gsl_vector* fw);
