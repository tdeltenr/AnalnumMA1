#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "find_norm.h"
#include "multigrid_solver.h"
#include "multigrid_struct.h"
#include "umfpk.h"

double Vec_Vec_prod(int n,double* a,double* b){
	double sum = 0.0;
	/* Compute the product a^t*b for 2 nx1 vectors*/
	
	for(int i =0; i<n;i++){
		sum += a[i]*b[i];
	}
	
	return sum;
}

void Mat_Vec_prod(int n ,int* ia, int* ja, double*a, double* p, double** Ap){
	// Compute the Ap product in CG method
	double sum;

	for(int i = 0; i < n; i++){
		sum = 0.0;
		for(int nz_id = ia[i]; nz_id < ia[i+1]; nz_id++){
				sum += a[nz_id]*p[ja[nz_id]];
		}
		(*Ap)[i] = sum;
	}
	return;
}

void Vector_sum(int n,double** vec1, double* vec2,double coef, double* vec3) {
    // Compute vec1 = vec2 + coef * vec3
        for (int i = 0; i < n; i++) {
        (*vec1)[i] = vec2[i] + coef * vec3[i];
    }
}

double CG_method(int it, int n,int m, int* ia, int* ja, double* a, double* b,double** u_m,int numLevels){
	// Compute the Conjugate Gradient method
	double r;
	double res = 1;
	double a_m,b_m,den;
	double* r_m = (double*)malloc(n*sizeof(double));
	double* invBr_m = (double*)malloc(n*sizeof(double));
	double* p_m = (double*)malloc(n*sizeof(double));
	double* d_m  = (double*)calloc(n,sizeof(double));
	double* Ad_m = (double*)malloc(n*sizeof(double));
	double* prev_r_m = (double*)malloc(n*sizeof(double));
	double* prev_invBr_m = (double*)malloc(n*sizeof(double));
	
	//Créer la grille pour le préconditionneur
	MultigridGrid** grids;
	grids = CreateMultiGridHierarchy(m,n,numLevels);
	
	// Factoriser la grille au coarsest level
	void *Numeric = NULL;	
	factorize_umfpack(grids[numLevels-1]->n,grids[numLevels-1]->ia,grids[numLevels-1]->ja,grids[numLevels-1]->a,grids[numLevels-1]->b,*u_m, &Numeric);
	
	//Créer le vecteur résidu
	residu_vector(n,ia,ja,a,b,*u_m,&r_m);
	
	
	// 1st iteration 
	b_m = 0;
	V_multigrid_preconditioning(grids[0],m,n,r_m,&invBr_m,Numeric); // calcule inBr_m
	Vector_sum(n,&d_m,invBr_m,b_m,d_m);
		
	Mat_Vec_prod(n,ia,ja,a,d_m,&Ad_m);
	a_m = Vec_Vec_prod(n,r_m,invBr_m)/Vec_Vec_prod(n,d_m,Ad_m);
	
	Vector_sum(n,u_m,*u_m,a_m,d_m);
	
	memcpy(prev_r_m,r_m,n*sizeof(double));
	memcpy(prev_invBr_m,invBr_m,n*sizeof(double));
	
	Vector_sum(n,&r_m,r_m,-a_m,Ad_m);

	r = residu(n,ia,ja,a,b,*u_m);
	printf("it : 0 = %e\n",r);
	
	for(int i =0; i<it-1 ;i++){
	//Préconditionner r_m
	V_multigrid_preconditioning(grids[0],m,n,r_m,&invBr_m,Numeric);
	b_m = Vec_Vec_prod(n,r_m,invBr_m)/Vec_Vec_prod(n,prev_r_m,prev_invBr_m);
	
	Vector_sum(n,&d_m,invBr_m,b_m,d_m);
	
	Mat_Vec_prod(n,ia,ja,a,d_m,&Ad_m);
	a_m = Vec_Vec_prod(n,r_m,invBr_m)/Vec_Vec_prod(n,d_m,Ad_m);
	
	Vector_sum(n,u_m,*u_m,a_m,d_m);
	
	
	Vector_sum(n,&r_m,r_m,-a_m,Ad_m);
	
	r = residu(n,ia,ja,a,b,*u_m);
	printf("it : %d = %e\n",i+1,r);
	}
	
	
	//Free les structures multigrilles
	for(int i = 0;i<numLevels;i++){
			freeMultigridHierarchy(grids[i]);
	}
	free(grids);
	free(Numeric);
	// Free les vecteurs
	free(r_m); free(invBr_m);free(p_m);free(d_m);free(Ad_m);free(prev_r_m);free(prev_invBr_m);
	
	return 0.0;

}




