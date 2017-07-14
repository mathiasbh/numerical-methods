typedef struct {int size; double* data;} nvector;

nvector* nvector_alloc       (int n);                           /* allocates memory for size-n vector */
void     nvector_free        (nvector* v);                      /* frees memory */
void     nvector_set         (nvector* v, int i, double value); /* v_{i} = value */
double   nvector_get         (nvector* v, int i);               /* returns v_{i} */
void     nvector_set_zero    (nvector* v);                      /* all elements = 0 */
double   nvector_dot_product (nvector* u, nvector* v);          /* dot-product */
void     nvector_add         (nvector* a, nvector* b);          /* a_i <- a_i + b_i */
void     nvector_sub         (nvector* a, nvector* b);          /* a_i <- a_i - b_i */
void     nvector_scale       (nvector* a, double x);            /* a_i <- x*a_i     */
