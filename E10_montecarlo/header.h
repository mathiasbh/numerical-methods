void samplingRand(int size, double* a, double* b, double *x);

void plainmc(int size, double f(double* x), double* a, double* b, int N, double* result, double* error);

void mctol(int size, double f(double* x), double* a, double* b, double tol, double* result, double* error);

double adaptOpen(double f(double), double a, double b,
		 double acc, double eps, double f2, double f3,
		 double nrec, double* err);

double intOpen(double f(double), double a, double b,
	       double acc, double eps, double* err);

double intOpen2d(int size, double f(double* x), double* a, double* b, double acc, double eps, double* error);
