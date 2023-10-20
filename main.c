#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "find_norm.h"
#include "plot.h"
#include "gs.h"
#include "multigrid_methods.h"
#include "multigrid_solver.h"
#include "multigrid.h"
#include "CG_Methods.h"


/* Fonction main */

int main(int argc, char *argv[])
{
  /* déclarer les variables */
  char solver_type;
  int m, it;
  int n, *ia, *ja; 
  double *a, *b,*res,r;
  double tc1, tc2, tw1, tw2, relative_error; /* mis à jour le 13/10/22 */
  m = 801;  // m doit etre de type 8n+1 avec n%2(k-1) == 0;
  it = 50;
  
  double* res_vector = malloc(it * sizeof(double)); 
  
  /* générér le problème */
  if (prob(m, &n, &ia, &ja, &a, &b))
     return 1;
   
  double* x_direct = calloc(n,sizeof(double)); 
       	
  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );

  //CG_method(it,n,m,ia,ja,a,b,&x_direct);
  
  //print_vector_double(n,x_direct);
  printf("2-GRID\n");
  two_grid_method(it,n,m,ia,ja,a,b,&res_vector);
  
  printf("Multi-Grid\n");
  V_multigrid(m,n,it,&x_direct,&res_vector);
  //plot_2D_graphs(it,res_vector);

   //allouer la mémoire pour le vecteur de solution 
  /*
  u_m = malloc(n * sizeof(double));
  if ( u_m == NULL) {
  	printf("\n ERREUR : pas de mémoire pour les vecteurs des solutions\n\n");
        return 1;
  }*/
  
  /*
  int k,m_c,n_c;
  k = (m-1)/8;
  m_c = k/2*8 +1; //Taille de la nouvelle grille
  n_c = (m_c-2)*(m_c-2) - (m_c-1)*(m_c-1)*9/64; // n pour la nouvelle grille
  
  
  double* u_c = malloc(n_c * sizeof(double));
  double* u_m = malloc(n * sizeof(double));
  double* r_c = malloc(n_c * sizeof(double));
  double* r_m = malloc(n * sizeof(double));
  
  for (int i = 0; i < n_c; i++) {
        u_c[i] = 1.0; // You can initialize to any desired value
    }
   
   for (int i = 0; i < n; i++) {
        r_m[i] = 1.0; // You can initialize to any desired value
    }
   
   fw_Restriction(m, n, r_m, &r_c);
   Prolongation(m,n,m_c,n_c,u_c,&u_m);

  
  
 /*
  for(int q = 0; q < 1; q++){
  r = backward_gauss_seidel(n,ia,ja,a,b,&x_direct);
  residu_vector(n,ia,ja,a,b,x_direct,&res);
  norm(n,res);
  fw_Restriction(m,n,res,&r_c);
  norm(n_c,r_c);
	}
  free(res);
  */
  
  /* Pour les prints : 
   * ia: int j =0; j< m*m - ((m-1)/2 + 1) - ((m-1)/4 +1) - (m-1)*(m-1)*9/64 %d
   * ja: int j=0; j < 5*n - 4*m - 4; j++ %d
   * a : int j=0; j < 5*n - 4*m - 4; j++ %f*/
  
  // Print les matrices si besoin
  //int k = m*m - 4*m -1 - (m-1)*(m-1)*4/64;
  /*
  FILE *f = NULL;
  f = fopen("ia.txt","w"); 
  for(int j=0; j< m*m - 4*m - (m-1)*(m-1)*4/64; j++){
	  fprintf(f,"%d \n",ia[j]);
  } 
  fclose(f);
  
   FILE *g = NULL;
  g = fopen("ja.txt","w"); 
  for(int j=0; j<5*k - 4*(m-2); j++){
	  fprintf(g,"%d \n",ja[j]);
  } 
  fclose(g);
  
   FILE *h = NULL;
  h = fopen("a.txt","w"); 
  for(int j=0;j<5*k - 4*(m-2); j++){
	  fprintf(h,"%f \n",a[j]);
  } 
  fclose(h);
  
  // SOLVEUR DIRECT
  /* résoudre et mesurer le temps de solution */
  
  /*
  tc1 = mytimer_cpu(); tw1 = mytimer_wall(); 
  if( solve_umfpack(n, ia, ja, a, b, x_direct) )
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall(); 
  printf("\nTemps de solution (CPU) [DIRECT]: %5.1f sec",tc2-tc1); 
  printf("\nTemps de solution (horloge) [DIRECT]: %5.1f sec \n",tw2-tw1);
   
  // Obtenir l'erreur relative
  relative_error = residu(n, ia,ja,a,b,x_direct); 
  printf("relative error [DIRECT] = %.16e \n", relative_error); 
  */

  
  //Plot le graphe
  //plot(m,x_direct);
 
  /* libérer la mémoire*/
  free(ia); free(ja); free(a); free(b); free(res_vector);
 
  return 0;
}

   
