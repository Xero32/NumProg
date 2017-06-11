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

void
setup_chebyshev_interpolationpoints(pinterpolation inter, double a, double b){
    int m = inter->m;
    double *x = inter->xi;
    assert(a<=b);
    assert(m);
    double z = (b-a) / 2;
    double w = (b+a) / 2;
    
    for(int i = 0; i < m+1; i++){
        x[i] = w + z * cos( (2*i - 1) * M_PI / 2 * m);
    }
}

void 
setup_aequidistant_interpolationpoints(pinterpolation inter, double a, double b){
    int m = inter->m;
    double *x = inter->xi;
    assert(b>=a);
    assert(m);
    double z = (double)(b-a)/m;
    printf("inside aequidistant fct\n");
    for(int i = 0; i < m+1; i++){
        x[i] = i * z - 1;
        printf("z: %f, x[i]: %f, i: %d\n",z,x[i],i);
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

double
n(int i, int j, double xx, pinterpolation inter){
    double *xi = inter->xi;
    double result = 1;
    
    for(int k = i; k < j; k++){
        result *= xx-xi[k];
    }
    return result;
}

void
newton_divided_differences(pinterpolation inter){
// A d = b  solve
// optimize!!
    int m = inter->m;
    double *x = inter->xi;
    double *y = inter->f;
    double *d = inter->d;
    double p[m+1];
    
    for(int i = 0; i < m; i++){
        for(int j = 0; j <= i; j++){
            p[j] = n(0,j,x[i], inter);
        }
        d[i] = y[i];
        for(int j = 0; j <= i; j++){
            d[i] -= p[j]*d[j];
        }
        assert(p[i]);
        d[i] /= p[i];
    }
}

/* Horner schema */

double
eval_interpolation_polynomial(const pinterpolation inter, double x){
    int m = inter->m;
    double *xi = inter->xi;
//     double *y = inter->f;
    double *d = inter->d;
    double p;
    
    p = d[m];
    
    for(int l = 1; l < m; l++){
        p = d[m-l] + (x - xi[m-l]) * p;
    }
    return p;
}
