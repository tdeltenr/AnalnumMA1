#include <stdio.h>
#include <stdlib.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "rho.h"
#include "find_norm.h"
#include "plot.h"
#include "flux.h"

/* Fonction main */

int main(int argc, char *argv[])
{
  /* déclarer les variables */
  char solver_type;
  int m;
  int n, *ia, *ja; 
  double *a, *b, *x_direct;
  double tc1, tc2, tw1, tw2, relative_error; /* mis à jour le 13/10/22 */
  m = 17;  // m doit etre de type 8n+1;
  
  /* générér le problème */
  if (prob(m, &n, &ia, &ja, &a, &b))
     return 1;
     	
  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );

  /* allouer la mémoire pour le vecteur de solution */

  x_direct = malloc(n * sizeof(double));
  if ( x_direct == NULL) {
  	printf("\n ERREUR : pas de mémoire pour les vecteurs des solutions\n\n");
        return 1;
  }
  
  /* Pour les prints : 
   * ia: int j =0; j< m*m - ((m-1)/2 + 1) - ((m-1)/4 +1) - (m-1)*(m-1)*9/64 %d
   * ja: int j=0; j < 5*n - 4*m - 4; j++ %d
   * a : int j=0; j < 5*n - 4*m - 4; j++ %f*/
  
  // Print les matrices si besoin
  int k = m*m - 4*m -1 - (m-1)*(m-1)*4/64;
  
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
  
  tc1 = mytimer_cpu(); tw1 = mytimer_wall(); 
  if( solve_umfpack(n, ia, ja, a, b, x_direct) )
     return 1;
  tc2 = mytimer_cpu(); tw2 = mytimer_wall(); 
  printf("\nTemps de solution (CPU) [DIRECT]: %5.1f sec",tc2-tc1); 
  printf("\nTemps de solution (horloge) [DIRECT]: %5.1f sec \n",tw2-tw1);
   
  // Obtenir l'erreur relative
  relative_error = residu(n, ia,ja,a,b,x_direct); 
  printf("relative error [DIRECT] = %.16e \n", relative_error); 


  
  //Plot le graphe
  //plot(m,x_direct);
 
  /* libérer la mémoire*/
  free(ia); free(ja); free(a); free(b); free(x_direct);
 
  return 0;
}

   
