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
  

/* Note:
    -'t' time value
    - 'delta' increment for the leapfrog method
    - 'data' used for m, c and the knowledge if a 
       wave on the boundary should start at the left
       interval edge or the right .
*/

  void
  step_leapfrog1d_wave(pcgridfunc1d u_old, pcgridfunc1d v_old,
  pgridfunc1d u_new, pgridfunc1d v_new, double t, double delta, void* data){
	
	/*
	- Was ist m und wofür ist Data genau gedacht?
	- void pointer = belibieger Datentyp?
	- Wo kriegen wir k her ? Analogon zu current ?
	- Warum gibt es d wenn es nur als n+2 definiert ist ?
	- Sind die Randwerte immer 0 ? Vergleich Aufgabenzettel(Veränderung der Randwerte)
	- Was bedeutet Left or Right wave aus Data ? Stimmt unsere Interpretation ?
	
	*/
	
	double* ux = u_old->x;
	double* uv = v_old->x;
	double* nx = u_new->x;
	double* nv = v_new->x;
	int n = u_old->n;
	double c = (double)data[1]; // depending on the position of c in data
	double lr = (double)data[0];
	double h = ux->g->h; 

	int k; // how to define k because we cannot access current in this function
	if(k==0){
		for(i = 0; i < n+2; i++){
			nx[i] = ux[i];
		}
		if(lr<0){left_boundary_gridfunc1d(ux, t)};
		else if(lr>0){right_boundary_gridfunc1d(ux, t)};
	} else if(k==1){
		for(i = 0; i < n+2; i++){
			nv[i] = uv[i];
		}
	} else if(k%2==0){
		for(i = 1; i < n+1; i++){
			nx[i] = ux[i] + 2.0*delta*uv[i];
		}
		if(lr<0){
			left_boundary_gridfunc1d(ux, t);
			nx[n+1]=ux[n+1];
		}
		else if(lr>0){
			right_boundary_gridfunc1d(ux, t);
			nx[0] = ux[0];	
		}	
	}
	else if(k%2==1){
		for(i = 1; i < n+1; i++){
			nv[i] = uv[i] + 2.0*delta*c*c*(ux[i-1]-2*ux[i]+ux[i+1])(h*h);
		}
		nv[0]=uv[0];
		nv[n+1]=uv[n+1];
	}
		
		
	
  }		
