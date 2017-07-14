/*  ---------------------------
   / n-dimensional vector program
   / NAME: Mathias Bojsen-Hansen
   / DATE: 10-02-2016 */

#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

typedef struct {int size; double* data;} nvector;

nvector* nvector_alloc (int n) { // returns nvector pointer
  nvector* v = (nvector*)malloc(sizeof(nvector)); // creates nvector named v of size nvector type casted as nvector
  v->size = n;
  v->data = (double*)malloc(n*sizeof(double)); // set data memory to n elements of double
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free (nvector* v) { 
  free(v->data); free(v); // free memory from v.data and v
}

void nvector_set (nvector* v, int i, double value) {
  assert(0<=i && i<(v->size));
  v->data[i] = value;
}

double nvector_get (nvector* v, int i) {
  assert(0<=i && i<(v->size));
  return v->data[i];
}

void nvector_set_zero (nvector* v) {
  for (int j=0; j<(v->size); j++) v->data[j] = 0;
}

double nvector_dot_product (nvector* u, nvector* v) {
  double result = 0;
  for (int j=0; j<(u->size)-1; j++) result+=(u->data[j])*(v->data[j]);
  return result;
}

void nvector_add (nvector* a, nvector* b) {
  for (int j=0; j<(a->size); j++) a->data[j]+=b->data[j];
}

void nvector_sub (nvector* a, nvector* b) {
  for (int j=0; j<(a->size); j++) a->data[j]-=b->data[j];
}

void nvector_scale (nvector* a, double x) {
  for (int j=0; j<(a->size); j++) a->data[j]*=x;
}
