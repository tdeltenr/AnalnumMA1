#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "find_norm.h"
#include "gs.h"
#include "multigrid_methods.h"
#include "plot.h"
#include "multigrid.h"
#include "multigrid_struct.h"
#include "umfpk.h"


double two_grid_method(int it,int n, int m,int* ia, int* ja, double* a, double* b, double** residual_vector){
	
	double r;
	int k,m_c,n_c;
	
	k = (m-1)/8;
	m_c = k/2*8 +1; //Taille de la nouvelle grille
	n_c = (m_c-2)*(m_c-2) - (m_c-1)*(m_c-1)*9/64;
	
	double* u_m = malloc(n * sizeof(double));
	double* r_m = malloc(n * sizeof(double));
	double* r_c = malloc(n_c * sizeof(double));
	double* u_c = malloc(n_c * sizeof(double));
	int *ia_c,*ja_c;
	double *a_c, *b_c;
	
	// Génération de la grille intermédiaire
	prob(m_c,&n_c,&ia_c,&ja_c,&a_c,&b_c);
	
	// Initialisation de u_m et u_c à 0
	for(int j=0;j<n;j++){
		u_m[j] = 0.0;
	}
	
	for(int l = 0; l<n_c;l++){
		u_c[l] = 0.0;
	}
	
	for(int i =0; i < it; i++){
		//Pre-smoothing
		r = backward_gauss_seidel(n,ia,ja,a,b,&u_m);
		
		
		// Création vecteur résidu 
		residu_vector(n,ia,ja,a,b,u_m,&r_m);
		fw_Restriction(m,n,m_c,n_c,r_m,&r_c);
			
		// Résolution du problème sur la coarsed grid
		solve_umfpack(n_c,ia_c,ja_c,a_c,r_c,u_c);
		
		// Prolongation
		Prolongation(m,n,m_c,n_c,u_c,&u_m);
		
		forward_gauss_seidel(n,ia,ja,a,b,&u_m);
		r = residu(n,ia,ja,a,b,u_m);
			
		//printf("%e\n",r);
		(*residual_vector)[i] = r;
		}
		
		free(u_c); free(r_c); free(r_m);
		return r;
}


void V_scheme(MultigridGrid* Grid, double** u_m){
	int* ia = Grid->ia;
	int* ja = Grid->ja;
	double* a = Grid->a;
	double* b = Grid->b;
	int n = Grid->n;

	MultigridGrid* Next_Grid = Grid->coarser;


	if (Grid->isCoarsestGrid){
	

	} else {
		// Appliquer le smoothing
		backward_gauss_seidel(n,ia,ja,a,b,u_m);
		
		//Envoyer le résidu au bon endroit dans la prochaine grille (dans b pour la dernière dans r pour les autres) 
		if (Grid->coarser != NULL && Next_Grid->isCoarsestGrid){
			
		residu_vector(n,ia,ja,a,b,(*u_m),&(Grid->residual));
		fw_Restriction(Grid->m,Grid->n,Next_Grid->m,Next_Grid->n,Grid->residual,&(Next_Grid->residual));
		
		
		}else{
		//residu_vector(n,ia,ja,a,b,(*u_m),&r_m);
		}
		
		// Appliquer le post-smoothing
		forward_gauss_seidel(n,ia,ja,a,b,u_m);
		
		
		
	}
	
	return;
	
}


int V_multigrid(int m,int n, int it){
	/* Initialiser la structure pour L = 4 */
	MultigridGrid* finestGrid, *intermediateGrid1, *intermediateGrid2, *coarsestGrid;
	
	double* u_m = (double*)calloc(n,sizeof(double));
	
	// G-1h
	double* r_h = (double*)malloc(sizeof(double) * n);
	
	finestGrid = fillMultigridGrid(m,n,0);
	finestGrid->residual = r_h;
	
	// G-2h
	int n2h,m2h;
	int k = (m-1)/8;
	m2h = k/2*8 +1;
	n2h = (m2h-2)*(m2h-2) - (m2h-1)*(m2h-1)*9/64;
	double* r_2h = (double*)malloc(sizeof(double) * n2h);

	intermediateGrid1 = fillMultigridGrid(m2h,n2h,1);
	intermediateGrid1->residual = r_2h;
		
	// G-4h 
	int n4h,m4h;
	k = (m2h-1)/8;
	m4h = k/2*8 +1;
	n4h = (m4h-2)*(m4h-2) - (m4h-1)*(m4h-1)*9/64;
	double* r_4h = (double*)malloc(sizeof(double) * n4h);
	 
	intermediateGrid2 = fillMultigridGrid(m4h,n4h,0);
	intermediateGrid2->residual = r_4h;

	// G-8h
	int n8h,m8h;
	k = (m4h-1)/8;
	m8h = k/2*8 +1;
	n8h = (m8h-2)*(m8h-2) - (m8h-1)*(m8h-1)*9/64;
	double* r_8h = (double*)malloc(sizeof(double) * n8h);
	
	coarsestGrid = fillMultigridGrid(m8h,n8h,1);
	coarsestGrid->residual = r_8h;
	
	// Relation entre les grilles
	finestGrid->finer = NULL;
	finestGrid->coarser = intermediateGrid1;
	
	intermediateGrid1->finer = finestGrid;
	intermediateGrid1->coarser = intermediateGrid2;
	
	intermediateGrid2->finer = intermediateGrid1;
	intermediateGrid2->coarser = coarsestGrid;
	
	coarsestGrid->finer = intermediateGrid2;
	coarsestGrid->coarser = NULL;
	int nnz = 5*n8h - 4*(m8h-2);
	
			
	// Commencer le V scheme
	for(int i = 0; i<10;i++){
	V_scheme(finestGrid,&u_m);
	}
	
	
	
	/* Free les Grids */
	freeMultigridHierarchy(finestGrid);
	freeMultigridHierarchy(intermediateGrid1);
	freeMultigridHierarchy(intermediateGrid2);
	freeMultigridHierarchy(coarsestGrid);
	
	return 0; 
}


