/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 5 - Simulation Loesung der Wellengleichung 	 */
/* ------------------------------------------------------------- */
/*	Autoren: Christina Boerst				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

  #include <assert.h>		
  #include "gridfunc1d.h"
  #include "leapfrog1d.h"		 
  #include <stdio.h>

/* Note:
    -'t' time value
    - 'delta' increment for the leapfrog method
    - 'data' used for m, c and the knowledge if a 
       wave on the boundary should start at the left
       interval edge or the right .
*/

  void  // change x AND v in one leapfrog call
  step_leapfrog1d_wave(pcgridfunc1d u_old, pcgridfunc1d v_old,
  pgridfunc1d u_new, pgridfunc1d v_new, double t, double delta, void* data){
	double *ox = u_old->x;
        double *ov = v_old->x;
        double *nx = u_new->x;
        double *nv = v_new->x;
        int n = u_old->d;
        double c = ((double*)data)[1]; // make sure c is saved in data[1]
        double lr = ((double*)data)[0];
        double h = u_old->g->h; 
        double k = 2.0 * c * c / h / h * delta;
        for(int i = 1; i < n-1; i++){
            nx[i] = ox[i] + 2.0 * delta * ov[i];

        }
        if(lr == 'l'){
            left_boundary_gridfunc1d(u_new,t);
        }else if(lr == 'r'){
            right_boundary_gridfunc1d(u_new,t);
        }
        for(int i = 1; i < n-1; i++){
            nv[i] = ov[i] +  k * (nx[i-1] - 2.0 * nx[i] + nx[i+1]);

        }
        
        

/*        FILE *f = fopen("checkfl7.txt","a");
        if(f){
            f = fopen("checkfl7.txt","a");
            for(int i = 0; i < n; i += 5){
                if(nx[i] && nv[i]) fprintf(f,"%f\t%f\n",nx[i],nv[i]);
            }
            fprintf(f,"\ntime: %f\tnewline\n\n",t);
            
        }
        if(f) fclose(f);   */     
}	
