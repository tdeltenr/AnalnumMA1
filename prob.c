#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double u(double x, double y){
	return sin(sqrt(x*x + y*y));
}

int prob(int m, int *n, int **ia, int **ja, double **a, double **b)
/*
   But
   ===
   Générer le système linéaire n x n 
                          
                             Au = b                                   

   qui correspond à la discrétisation sur une grille cartésienne 
   régulière m x m de l'équation de Poisson à deux dimensions
              
            d    d        d    d
         - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
           dx   dx       dy   dy

  avec les conditions aux limites de Dirichlet
         
         u = 20     sur (0,y)                 , avec 0.5 <=  y  <= 1.5
         du/dn = 1  sur (1,y), (x,0) et (x,1) , avec 0 <= x,y <= 1 .
         u = 0      sur (x,0)                 , avec 1.5 <= x <= 3.5

  La numérotation des inconnues est lexicographique, la direction x étant 
  parcourue avant celle de y. La matrice A est retournée dans le format 
  CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.

  Arguments
  =========
  m  (input)  - nombre de points par direction dans la grille
                (les valeurs m inférieures à 2 ne sont pas valides) 
  n  (output) - pointeur vers le nombre d'inconnus dans le système
  ia (output) - pointeur vers le tableau 'ia' de la matrice A
  ja (output) - pointeur vers le tableau 'ja' de la matrice A
  a  (output) - pointeur vers le tableau 'a' de la matrice A
  b  (output) - pointeur vers le tableau 'b'

*/
{
    int  nnz, ix, iy, ind, next_ind, ub_ix;
    double h,invh2;

    if(((m-1)%8) != 0) { // Vérifier que le pas est une valeur valide !!! 
        printf("\n ERREUR : m = %d n'est pas une valeur valide\n\n",m);
        return 1;
    }
    h = 4/((double)m-1); //Cast double pour pas avoir  0
    invh2 = (m-1)*(m-1)/16; /* h^-2 pour L=4 */
    *n  = (m-2)*(m-2) - (m-1)*(m-1)*9/64;//m*m - (m-1)*(m-1)*4/64; /* nombre d'inconnues */
    nnz = 5*(*n) - 4*(m-2); /* nombre d'éléments non nuls 5n pour les pts intérieurs - 4m pour les murs - 4 pour les bords avec les fenêtres */
    
    /* allocation des tableaux */
    *ia  = (int*)malloc((*n + 1) * sizeof(int));
    *ja  = (int*)malloc(nnz * sizeof(int));
    *a   = (double*)malloc(nnz * sizeof(double));
    *b   = (double*)malloc(*n*sizeof(double));
    /* allocation réussite? */

    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }

    /* partie principale : remplissage de la matrice */
	
    ind = 0; /* au cas où m<=1 */
    next_ind = 0;
    nnz = 0; 
    int j = 0;
    for (iy =1 ; iy < m-1; iy++) { // y va bien de 1 à m-1 vu que les murs sont en "Dirichlet" 
			if(iy >= (m-1)*5/8){ // Définir l'upper bound sur ix
				ub_ix = (m-1)*5/8;
			} else{
				ub_ix = m-1;
			}
		for (ix =1;ix < ub_ix; ix++){
				/* Trouver l'indice du point et l'incrémenter pour la suite */
				ind = next_ind;
				next_ind++;
				/* marquer le début de la ligne suivante dans le tableau 'ia' */
				(*ia)[ind] = nnz;
				
				/* (v) remplissage de la ligne : voisin sud */
				if (iy == 1) { /* Condition Dirichlet u du bord sud */         
					(*b)[ind] += u(h*ix,h*(iy-1))*invh2;
				} else { /* reste du domaine */
					(*a)[nnz] = -invh2; 
					if (iy <= (m-1)*5/8){
					(*ja)[nnz] = ind - m + 2 ; // Voisin sud au bon endroit 
					}else {
						(*ja)[nnz] = ind - (m-1)*5/8 + 1 ; //Voisin sud plus proche en terme d'indice 
					}
					nnz++;
				} 

				/* remplissage de la ligne : voisin ouest */
				if(ix == 1) { /* Dirichlet sur u à l'ouest */
					(*b)[ind] += u(h*(ix-1),h*iy)*invh2;
					
				} else { /* milieu du domaine */
					(*a)[nnz] = -invh2; 
					(*ja)[nnz] = ind - 1;
					nnz++;
				}
			
     
				/* remplissage de la ligne : élément diagonal */
				(*a)[nnz] = 4.0*invh2;
				(*ja)[nnz] = ind;
				nnz++;
	

				/* remplissage de la ligne : voisin est */
				if ((ix == m-2 && iy < (m-1)*5/8 )|| (iy >= (m-1)*5/8 && ix == (m-1)*5/8-1)){ /* Condition Dirichlet u à l'est */
					(*b)[ind] += u(h*(ix+1),h*iy)*invh2;

				} else if ((iy <= (m-1)*5/8 && ix < m-1) || (iy > (m-1)*5/8 && ix <= ((m-1)*5/8)-1) )  { /* milieu du domaine */
					(*a)[nnz] = -invh2;
					(*ja)[nnz] = ind + 1;
					nnz++;
				} 

				/* remplissage de la ligne : voisin nord */
				if (iy == m-2 && ix < (m-1)*5/8 || iy == (m-1)*5/8-1 && ix >= (m-1)*5/8) { /* condition Dirichlet u du bord nord */         
					(*b)[ind] += u(h*ix,h*(iy+1))*invh2;
	
				} else if ((ix <= (m-1)*5/8 && iy < m-1) || (ix > (m-1)*5/8 && iy < (m-1)*5/8) ) { /* reste des points */
					(*a)[nnz] = -invh2;
					if (iy < (m-1)*5/8){
					(*ja)[nnz] = ind + m - 2;
				   } else{
					(*ja)[nnz] = ind + (m-1)*5/8 - 1;
					}
					nnz++;
            }
         }
      }
	
    /* dernier élément du tableau 'ia' */
    (*ia)[ind + 1] = nnz;
    
    /* retour habituel de fonction */
    return 0;
}
