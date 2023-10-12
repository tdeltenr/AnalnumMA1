#include <stdlib.h>
#include <stdio.h>

int plot(int m, double *x) {
	/* Affiche le graphe via gnuplot */
	int ind = 0,ub_ix,ix,iy;
	double h;
	double border_value = 0.0;
    
    h = 4.0/(m-1);
	FILE *f = NULL;
	f = fopen("data.txt","w"); 
	
	 for (iy =0 ; iy < m-1; iy++) { // y va bien de 1 à m-1 vu que les murs sont en "Dirichlet" 
			if(iy >= (m-1)*5/8){ // Définir l'upper bound sur ix
				ub_ix = (m-1)*5/8;
			} else{
				ub_ix = m-1;
			}
		for (ix =0;ix < ub_ix; ix++){
			fprintf(f, "%f %f ",ix*h,iy*h);
			//Gèrer les cas où on doit afficher une autre valeur que celle venant de x;
			if(iy == 0){
				fprintf(f,"%f",border_value);	
			}else if(ix == 0 || iy ==0){
				fprintf(f,"%f",border_value);
			}else if(ix == ub_ix){
				fprintf(f,"%f",border_value);
			}else if((iy == m-1 && ix <= (m-1)*5/8)||(iy == (m-1)*5/8 && ix > (m-1)*5/8)){
				fprintf(f,"%f",border_value);
			
			}else {
				fprintf(f,"%f",x[ind]);
				ind++;
			}
			fprintf(f,"\n");
		}
		fprintf(f,"\n");
	}
	fclose(f);
	system("gnuplot -p gnuplot.gp\\  ");
	return 1;
}

void plot_2D_graphs(int it, double* vector){
	FILE *f = NULL;
	f = fopen("data.txt","w");
	
	for(int i = 0; i < it; i++){
		fprintf(f, "%d %e \n ",i+1,vector[i]);
	}
	fclose(f);
	system("gnuplot -p gnuplot.gp\\  ");
	
	return;
}
	
