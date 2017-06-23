/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 8 - Approximation von Integralen Quadratur  	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker			 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

  #include <stdio.h>	// TODO remove
  #include <stdlib.h>
  #include <assert.h>		
  #include <math.h>
  #include "quadrature.h"

#ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif		
  
#define EPS 1.0e-10 
 /* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pquadrature
new_quadrature(int m){

  pquadrature quad;

  quad = (pquadrature) calloc(1, sizeof(quadrature));
  quad->xq = (double *) calloc(m+1, sizeof(double));
  quad->w = (double *) calloc(m+1, sizeof(double));
  quad->m = m;
  
  return quad;
}

void
del_quadrature(pquadrature quad){

  free(quad->xq);
  free(quad->w);
  free(quad);

}

pquadrature
setup_midpointrule(void){
    pquadrature q = new_quadrature(0);
    double *xq = q->xq;
    double *w = q->w;
    
    
    xq[0] = 0.0;
    w[0] = 2.0;
    return q;
}

pquadrature
setup_trapezoidalrule(void){
    pquadrature q = new_quadrature(1);
    double *xq = q->xq;
    double *w = q->w;
    
    xq[0] = -1.0;
    xq[1] = 1.0;
    w[0] = w[1] = 1.0;
    return q;
}

void
map_quadrature_points(pquadrature quad, double a, double b, double *x){
    double *xq = quad->xq;
    int m = quad->m;
    double d = (b-a) * 0.5;
    double h = (b+a) * 0.5;

    for(int i = 0; i < m+1; i++){
        x[i] = h + d * xq[i];
}



double
eval_quadrature(pquadrature quad, double a, double b, function f, void *data){
    int m = quad->m;
    double *x;
    x = calloc(m+1, sizeof(double));
    double *w = quad->w;
    
    double d = (b-a) * 0.5;
    double A = 0.0;
    
    
    map_quadrature_points(quad,a,b,x);
    
    for(int i = 0; i < m+1; i++){
        A += w[i] * f(x[i],data);
    }
    free(x);
    return  A * d;
}

double
eval_composite_quadrature(pquadrature quad, double a, double b, int n, function f, void *data){
    assert(n);
    double l = b-a;
    double u = 1.0/(double)n;
    double result = 0.0;
    for(int i = 0; i < n; i++){            
        result += eval_quadrature(quad, (double)i*l*u+a, (double)(i+1)*l*u+a, f, data);
    }
    return result;
}
