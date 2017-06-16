/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 7 - Approximation von Funktionen Interpolation  	 */
/* ------------------------------------------------------------- */
/*	Autoren: 						 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

  #include <stdlib.h>
  #include <assert.h>		
  #include <math.h>
  #include "interpolation.h"

#ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif		
  
 /* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pinterpolation
new_interpolation(int m){

  int i;
  pinterpolation inter;

  inter = (pinterpolation) calloc(1, sizeof(interpolation));
  inter->xi = (double *) calloc(m+1, sizeof(double));
  inter->f = (double *) calloc(m+1, sizeof(double));
  inter->d = (double *) calloc(m+1, sizeof(double));
  inter->m = m;
  
  for(i = 0; i < m+1; i++){
    inter->d[i] = 0.0;
  }

  return inter;
}

void
del_interpolation(pinterpolation inter){

  free(inter->xi);
  free(inter->f);
  free(inter->d);
  free(inter);

}

/* Chebyshev points [a,b] */

void
setup_chebyshev_interpolationpoints(pinterpolation inter, double a, double b){

	int m = inter->m;
	double* xi = inter->xi;
	assert(a<=b);
	double rad = (b-a)/2;
	double mid = (a+b)/2;
	
	for(int i=0 ; i<m+1;i++){
		xi[i] = mid + rad*cos((double)(2*i+1)/(2*m+2)*M_PI);
	}

}

void 
setup_aequidistant_interpolationpoints(pinterpolation inter, double a, double b){

	int m = inter->m;
	double* xi = inter->xi;
	assert(a<=b);
	double rad = (b-a)/m;
	
	for(int i=0 ; i<m+1;i++){
		xi[i] = b - rad*i;
	}
	
}

void
eval_interpolated_values(pinterpolation inter, function f, void *data){
	
	int m = inter->m;
	double* xi = inter->xi;
	double* yi = inter->f;
	
	for(int i=0;i<m+1;i++){
		yi[i] = f(xi[i],data);	
	}
	
}

double
n(int i,int j,double x, const pinterpolation inter){
	
	double* xi = inter->xi;	
	double n=1;
	
	for(int k=i;k<j;k++){
		n *= x-xi[k];
	}
	
	return n;
}

/* Evaluate Newton divided differences */

void
newton_divided_differences(pinterpolation inter){
	
	int m = inter->m;
	double* xi = inter->xi;
	double* yi = inter->f;
	double* di= inter->d;
	
	di[0]=yi[0];
	
	for(int k=1; k<m+1; k++){
		double n=1;
		di[k]=0.0;
		for(int i=0; i<k; i++){
			di[k]-=n*di[i];
			n*=xi[k]-xi[i];
		}
		di[k]+=yi[k];
		di[k]/=n;
	}
	
}

/* Horner schema */

double
eval_interpolation_polynomial(const pinterpolation inter, double x){

	int m = inter->m;
	double* xi = inter->xi;
	double* di= inter->d;
	double p = di[m];
	
	for(int l=1;l<=m;l++){
		p *= x-xi[m-l];
		p += di[m-l];
	}
	
	return p;
}

