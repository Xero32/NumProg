/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/*   	Serie 10 - Schwachbesetzte Matrizen (Richardson)	 */
/* ------------------------------------------------------------- */
/*	Autoren: 		 				 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

  #include <stdio.h>		
  #include <stdlib.h> 
  #include <math.h>
  #include "crs.h"

  #ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif

/* right hand side */
static double 
sinus(double x, void *data){
	
  double *dat = (double *) data;
  double fx;	
  double mu = dat[0];
	  
  fx = sin(M_PI * mu * x);	
  
  return fx;
}


int 
main(void){
  
  pcrsmatrix crs;
  unsigned int m;
  pvector x, b;
  unsigned int i;
  double data[1];
  double theta, h, lambda;

  m = 5;					/* problem dimension */
  data[0] = 1.0;
  theta = 0.01;
  h = 1.0/(m+1);
  
 
  printf("Solve poisson equation with Richardson iteration\n");
   printf("-----------------------------------------\n");
 
  /* set up problem */  
  crs = setup_poisson(m);  	/* matrix */
  b = new_zero_vector(m);	/* right hand side */
  x = new_zero_vector(m);	/* solution */
  set_righthandside(m, sinus, data, b);
  
  printf("Right hand side \n");
    for(i = 0; i < m; i++){
    printf("%f\n", b->x[i]);
  }
  /* solve with Richardson iteration */
  richardson_iteration(crs, x, theta, b, 1.e-7);
  printf("-----------------------------------------\n");
  printf("Solution\n");
  	lambda = sinus(h *0.5, data);		/* sin(k*h/2)*/
    lambda *= lambda;					/* sin^{2} */
 	lambda *= 4.0/h/h;					/* 4.0/(h*h) * sin */
  for(i = 0; i < m ; i++){
    printf("Richardson %f \t exakt %f \t difference %e \n", x->x[i],  b->x[i]/ lambda, fabs(x->x[i] - b->x[i]/ lambda));
  }
 
  /* Cleaning up */
  del_crsmatrix(crs);
  del_vector(x);
  del_vector(b);

  return EXIT_SUCCESS;
}
