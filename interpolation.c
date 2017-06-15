/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 7 - Approximation von Funktionen Interpolation  	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker         	 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

  #include <stdlib.h>
  #include <assert.h>		
  #include <math.h>
  #include "interpolation.h"
#include <stdio.h>
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
// m-te Ordnung bedeutet m+1 Interpolationspunkte! Korrektur!
void
setup_chebyshev_interpolationpoints(pinterpolation inter, double a, double b){
    int m = inter->m;
    double *x = inter->xi;
    if(a > b){
        double c = a;
        a = b;
        b = c;
    }
    assert(m);
    double z = (b-a) * 0.5;
    double w = (b+a) * 0.5;
    // stattdessen lieber Multiplikation
    for(int i = 0; i < m+1; i++){
        x[m-i] = w + z * cos( (2.0 * (double)i + 1.0) * M_PI * 0.5 / ((double)m + 1.0) );
    }
}

void 
setup_aequidistant_interpolationpoints(pinterpolation inter, double a, double b){
    int m = inter->m;
    double *x = inter->xi;
    if(a > b){
        double c = a;
        a = b;
        b = c;
    }
    assert(m);
    double z = (double)(b-a)/(m);
    for(int i = 0; i < m+1; i++){
        x[i] = i * z + a;
    }
}

void
eval_interpolated_values(pinterpolation inter, function f, void *data){
    int m = inter->m;
    double *x = inter->xi;
    double *y = inter->f;
    
    for(int i = 0; i < m+1; i++){
        y[i] = f(x[i], data);
    }
}

/* Evaluate Newton divided differences */
// d[m] falsch
// einfache methode, um d[i] zu checken ?
void
newton_divided_differences(pinterpolation inter){
    int m = inter->m;
    double *x = inter->xi;
    double *y = inter->f;
    double *d = inter->d;
    d[0] = y[0];
    for(int u = 0; u <= m; u++){
        d[u] = y[u];
    }
    for(int k = 1; k <= m; k++){
        for(int j = m+1; j-->k;){ // int j = m+1
            int i = j-k;
            assert(x[j] - x[i]);
            d[j] = (d[j] - d[j-1]) / (x[j] - x[i]);
        }
    }
}

/* Horner schema */

double
eval_interpolation_polynomial(const pinterpolation inter, double x){
    int m = inter->m;
    double *xi = inter->xi;
    double *d = inter->d;
    double p = d[m];
    
    for(int l = 1; l <= m; l++){
        p *= x-xi[m-l];
        p += d[m-l];
    }

    return p;
}