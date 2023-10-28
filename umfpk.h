/* prototype */
int solve_umfpack(int n, int *ia, int *ja, double *a, 
                  double *b, double *x,double* Numeric);
                  
int factorize_umfpack(int n, int *ia, int *ja, double *a, double *b, double *x,void **NumericOut);
