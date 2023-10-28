#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "multigrid.h"
#include "find_norm.h"
#include "prob.h"

MultigridGrid* fillMultigridGrid(int m, int n, int isCoarsestGrid) {
    MultigridGrid* grid = (MultigridGrid*)malloc(sizeof(MultigridGrid));
	
    // Declare pointers for ia, ja, a, and b
    int* ia;
    int* ja;
    double* a;
    double* b;
	int nnz = 5*n - 4*(m-2);
    // Call prob to fill the arrays
    prob(m, &n, &ia, &ja, &a, &b);
	
	
    // Allocate memory for new arrays and copy data
    grid->m = m;
    grid->n = n;
    grid->isCoarsestGrid = isCoarsestGrid;
    grid->ia = (int*)malloc(sizeof(int) * (n+1));
    grid->ja = (int*)malloc(sizeof(int) * nnz);
    grid->a = (double*)malloc(sizeof(double) * nnz);
    grid->b = (double*)malloc(sizeof(double) * n);
    grid->residual = (double*)malloc(sizeof(double) * n);
	
	
    // Copy data from the returned arrays
    memcpy(grid->ia, ia, sizeof(int) * (n+1));
    memcpy(grid->ja, ja, sizeof(int) * nnz);
    memcpy(grid->a, a, sizeof(double) * nnz);
    memcpy(grid->b, b, sizeof(double) * n);
       
    // Free the arrays returned by prob
    free(ia);
    free(ja);
    free(a);
    free(b);
  
    return grid;
}


MultigridGrid** CreateMultiGridHierarchy(int m, int n, int numLevels){
	MultigridGrid** grids = (MultigridGrid**)malloc(sizeof(MultigridGrid*) * numLevels);
	
	// Rempli le niveau le plus fin de la grille
	int mFinest = m;
	int nFinest = n;
	double* rFinest = (double*)calloc(n,sizeof(double));
	grids[0] = fillMultigridGrid(mFinest, nFinest, 0);
	grids[0]->residual = rFinest;
	grids[0]->finer = NULL;

	// Créez les niveaux de grille restants
	for (int level = 1; level < numLevels; level++) {
		int k = (mFinest - 1) / 8;
		mFinest = k / 2 * 8 + 1;
		nFinest = (mFinest - 2) * (mFinest - 2) - (mFinest - 1) * (mFinest - 1) * 9 / 64;
		double* r = (double*)calloc(nFinest,sizeof(double));

		grids[level] = fillMultigridGrid(mFinest, nFinest, 0);
		grids[level]->residual = r;
		grids[level]->finer = grids[level-1];
		grids[level-1]->coarser = grids[level];
	}

	// Configurez la relation pour le niveau le plus grossier
	grids[numLevels - 1]->coarser = NULL;
	grids[numLevels - 1]->isCoarsestGrid = 1;
	
	return grids;
}

void freeMultigridHierarchy(MultigridGrid* grid) {
    // Free les éléments dans la Multigrid
    free(grid->ia);
    free(grid->ja);
    free(grid->a);
    free(grid->b);
    free(grid->residual);

    // free la structure en elle même 
    free(grid);
}
