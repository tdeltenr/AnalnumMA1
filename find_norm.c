#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double residu(int n,int *ia, int *ja, double *a, double *b, double *x){
	/* Renvoie la norme du résidu*/
	double Ax, r, normr, normb;
	normr = 0.0;
	normb = 0.0;
	
	for(int i = 0; i < n; i++){
		Ax = 0.0;
		for(int nz_id = ia[i]; nz_id < ia[i+1]; nz_id++){
			Ax += a[nz_id]*x[ja[nz_id]];
		}
		r = Ax - b[i];
		normr += r*r;
		normb += b[i]*b[i];
	}
	return sqrt(normr/normb);
}

void residu_vector(int n,int *ia, int *ja, double *a, double *b, double *x, double **res){
	double Ax;
	/*Remplis le vecteur résidu*/
	for(int i = 0; i < n; i++){
		Ax = 0.0;
		for(int nz_id = ia[i]; nz_id < ia[i+1]; nz_id++){
			Ax += a[nz_id]*x[ja[nz_id]];
		}
		(*res)[i] = b[i]-Ax;
	}	
	
	return;
}

double norm(int n, double *vector){
	double norme = 0.0;
	for(int i = 0; i<n; i++){
		norme += vector[i]*vector[i];
	}
	
	return sqrt(norme);
}

void print_vector_double(int n, double* vector){
	printf("[");
	for(int i=0;i<n;i++){
			if(i !=0 && i%10 == 0){
				printf("%f\n",vector[i]);
			} else {
				printf("%f,",vector[i]);
	}
	}
	printf("]");
	printf("\n");
	return;
}

void print_vector_int(int n, int* vector){
	printf("[");
	for(int i=0;i<n;i++){
			if(i !=0 && i%10 == 0){
				printf("%d\n",vector[i]);
			} else {
				printf("%d,",vector[i]);
	}
	}
	printf("]");
	printf("\n");
	return;
}
	
