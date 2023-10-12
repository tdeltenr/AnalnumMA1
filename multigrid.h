#include <stdio.h>
#include <stdlib.h>

/* Creation de la structure pour la multigrille */

typedef struct {
	int m; 
	int n;
	int isCoarsestGrid; // 1 for yes 0 for no
	
	//Pointeur vers les grilles au dessus et en dessous
	struct MultigridGrid* finer;
	struct MultigridGrid* coarser;
	
	int* ia;          
    int* ja;       
    double* a; 
    double* b; 
    double* residual;
	
} MultigridGrid;





