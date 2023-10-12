#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "multigrid.h"


void i_Restriction(int m, int n,int m_c, int n_c, double *r_m,double **r_c){
	// Reçoit le résidu sur G1h et renvoie un résidu sur la G2h avec une simple injection 
	int ix,iy,ub_ix,col_skip_j;
		
	int j = (m-2) +1; // Indice du premier élément de la coarse grid dans l'ancienne grille 
	int i = 0;
	
	for (iy =1 ; iy < m_c-1; iy++) { // On itère directement sur la coarsed grid 
			if(iy >= (m_c-1)*5/8){ 
				ub_ix = (m_c-1)*5/8;
				col_skip_j = (m-1)*5/8;
			}else{
				ub_ix = m_c-1;
				col_skip_j = m-1;
			}
				for (ix =1;ix < ub_ix; ix++){
					(*r_c)[i] = r_m[j];
					//printf("%f\n",(*r_c)[i]);
					j += 2;
					i += 1;
				}
			j += col_skip_j;
		}

	return;
}

void fw_Restriction(int m, int n,int m_c, int n_c, double *r_m, double **r_c){
	// Reçoit le résidu sur G1h et renvoie un résidu sur la G2h avec une simple injection 
	int ix,iy,ub_ix,col_skip_j,m_up,m_down;
	
	int j = (m-2) +1; // Indice du premier élément de la coarse grid dans l'ancienne grille 
	int i = 0;
	
	for (iy =1 ; iy < m_c-1; iy++) { // On itère directement sur la coarsed grid 
			if(iy >= (m_c-1)*5/8){ 
			ub_ix = (m_c-1)*5/8;
			col_skip_j = (m-1)*5/8;
			m_up = (m-1)*5/8-1; 
			if (iy == (m_c-1)*5/8){
				m_down = m-2;
			}else{
				m_down = (m-1)*5/8-1;
			}
		}else{
			ub_ix = m_c-1;
			col_skip_j = m-1;
			m_up = m-2;
			m_down = m-2;
			} 
				for (ix =1;ix < ub_ix; ix++){
					(*r_c)[i] = (r_m[j] + (r_m[j-1]+r_m[j+1]+r_m[j-m_down]+r_m[j+m_up])*1/2 + (r_m[j+m_up+1]+r_m[j+m_up-1]+r_m[j-m_down-1]+r_m[j-m_down+1])*1/4)*1/4;
							
					j += 2;
					i += 1;
				}	
				
			j += col_skip_j;
		}
			
	return;
}

void Prolongation(int m, int n,int m_c, int n_c, double *u_c, double **u_m){
	// Reçoit le vecteur solution de la G2h et le projette sur la G1h
	
	int ub_ix,col_skip_j,m_up,m_down,ix,iy,i = 0;
	int j = (m-2) +1;
	double uci;
	
	for (iy =1 ; iy < m_c-1; iy++) { // On itère directement sur la coarsed grid 
		if(iy >= (m_c-1)*5/8){ 
			ub_ix = (m_c-1)*5/8;
			col_skip_j = (m-1)*5/8;
			m_up = (m-1)*5/8-1; 
			if (iy == (m_c-1)*5/8){
				m_down = m-2;
			}else{
				m_down = (m-1)*5/8-1;
			}
		}else{
			ub_ix = m_c-1;
			col_skip_j = m-1;
			m_up = m-2;
			m_down = m-2;
			} 
	for (ix =1;ix < ub_ix; ix++){
			// Pour chaque point de la G2h il y a 8 points de la G1h qui l'entourent + le point lui même 
			uci = u_c[i];
			(*u_m)[j] += uci;
			(*u_m)[j+1] += uci*1/2;
			(*u_m)[j-1] += uci*1/2;
			(*u_m)[j+m_up] += uci*1/2;
			(*u_m)[j-m_down] += uci*1/2;
			(*u_m)[j+m_up-1] += uci*1/4;
			(*u_m)[j+m_up+1] += uci*1/4;
			(*u_m)[j-m_down-1] += uci*1/4;
			(*u_m)[j-m_down+1] += uci*1/4;
			 
			j += 2;
			i += 1;
	}
	j += col_skip_j;

	}	
	return;
}




