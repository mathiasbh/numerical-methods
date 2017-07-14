/* Perform binary search */

int binary_search(int n, double *x, double z) {
  int IndexInitial = 0; 
  int IndexFinal = n-1;

  while(IndexFinal - IndexInitial > 1) {
    int IndexMid = (IndexInitial + IndexFinal)/2;

    if(z > x[IndexMid]) IndexInitial = IndexMid; 
    else IndexFinal = IndexMid;
  }

  return IndexInitial;
}
