#include <stdlib.h>
#define RND ((double)rand()/RAND_MAX)
/* Distribute points uniformly through integration region */
void samplingRand(int size, double* a, double* b, double *x) {
  for(int i=0; i<size; i++) x[i] = a[i]+RND*(b[i]-a[i]);
}
