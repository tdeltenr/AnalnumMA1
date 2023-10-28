double Vec_Vec_prod(int n,double* a,double* b);

void Mat_Vec_prod(int n ,int* ia, int* ja, double*a, double* p, double** Ap);

void Vector_sum(int n,double** vec1, double* vec2,double coef, double* vec3);

double CG_method(int it, int n,int m, int* ia, int* ja, double* a, double* b,double** u_m,int numLevels);
