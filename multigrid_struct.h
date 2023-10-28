MultigridGrid* fillMultigridGrid(int m, int n, int isCoarsestGrid);

MultigridGrid** CreateMultiGridHierarchy(int m, int n, int numLevels);

void freeMultigridHierarchy(MultigridGrid* grid);

