/* Prototype */
double residu(int n, int *ia, int *ja, double *a, double *b, double *x); 

void residu_vector(int n,int *ia, int *ja, double *a, double *b, double *x, double **res);

double norm(int n, double *vector);

void print_vector_double(int n, double* vector);

void print_vector_int(int n,int* vector);
