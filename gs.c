#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "find_norm.h"

double forward_gauss_seidel(int n,int *ia, int *ja, double *a, double *b, double **x){
	
	// Reçois la matrice A et le vecteur x et utilise Gauss-Seidel
	double a_ii,r;
	
	// Une itération de forward Gauss Seidel
	for(int i = 0; i < n; i++){
		for(int nz_id = ia[i]; nz_id < ia[i+1]; nz_id++){
			if (ja[nz_id] != i){
				(*x)[i] -=   a[nz_id]*(*x)[ja[nz_id]];
			} else {
				a_ii = a[nz_id];
			}
		}
		(*x)[i] += b[i];
		(*x)[i] *= (double)1/a_ii;
	}
	
	r = residu(n,ia,ja,a,b,*x);
	return r;
}
