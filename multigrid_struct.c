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
