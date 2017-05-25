/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 5 - Simulation Loesung der Wellengleichung 	 */
/* ------------------------------------------------------------- */
/*	Autoren: Christina Boerst				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include "gridfunc1d.h"
  
  #ifndef M_PI
/* Copied from math.h*/
#define M_PI		3.14159265358979323846
#endif

/*-------------------------------------------------
*	 Create and destroy a grid 
*--------------------------------------------------*/


  pgrid1d
  new_grid1d(unsigned int n){

    pgrid1d g;
   
    g = calloc(1, sizeof(grid1d));
    g->n = n;
    g->h = 1.0/(n+1.0);
  
    return g;
  }

  void
  del_grid1d(pgrid1d g){

    free(g);

  }

/*------------------------------------------------
*      Create and destroy a grid function
*-------------------------------------------------*/

pgridfunc1d
new_gridfunc1d(pcgrid1d g){
	pgridfunc1d u_h;
	unsigned int n = g->n;
	
	u_h = calloc(1, sizeof(gridfunc1d));
	u_h->g = g;
	u_h->d = n+2;
	u_h->x = calloc(n+2, sizeof(double) );
	
	return u_h;
}

  void
  del_gridfunc1d(pgridfunc1d u_h){

    free(u_h->x);
    free(u_h);

  }

/*------------------------------------------------
*	    Utility functions
*-------------------------------------------------*/
void
zero_gridfunc1d(pgridfunc1d u_h){
  
    unsigned int i;
    unsigned int d = u_h->d;
    double *ux = u_h->x;


   for(i = 0; i < d; i++){
	ux[i] = 0.0;
  }
	
  }

void
left_boundary_gridfunc1d(pgridfunc1d u_h, double t){

	double *ux = u_h->x;

		if(t - 0.25 < 0 && 0 < t){
			ux[0] = sin(M_PI * (- t)*8.0);
		}
		else{
			ux[0] = 0.0;
		}

}

void
right_boundary_gridfunc1d(pgridfunc1d u_h, double t){

 	pcgrid1d g = u_h->g;
	double *ux = u_h->x;
	unsigned int n = g->n;

		if(t - 0.25 < 0 && 0 < t){
			ux[n+1] = sin(M_PI * (- t)*8.0);
		}
		else{
			ux[n+1] = 0.0;
		}

}

