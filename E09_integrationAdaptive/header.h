
double adaptOpen(double f(double), double a, double b,
		 double acc, double eps, double f2, double f3,
		 double nrec, double* err);

double intOpen(double f(double), double a, double b,
	       double acc, double eps, double* err);

double adaptClosed(double f(double), double a, double b,
		   double acc, double eps, double f1, double f3,
		   double nrec, double* err);

double intClosed(double f(double), double a, double b,
		 double acc, double eps, double* err);

double clensshawcurtis(double f(double), double a, double b,
		       double acc, double eps, double *err);
