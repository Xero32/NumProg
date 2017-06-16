
  #ifndef INTERPOLATION_H
  #define INTERPOLATION_H
  
  #include "basic.h"

  /* prototyp function*/
  typedef double (*function)(double x, void *data);

  typedef struct _interpolation interpolation;

  typedef interpolation *pinterpolation;
  
 
  struct _interpolation {
    int m;			/* Interpolation order */

    double *xi;			/* Interpolation points */
    double *f;			/* Interpolated values */
    double *d;			/* Newton divided differences */
  };

  
 /* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pinterpolation
new_interpolation(int m);

void
del_interpolation(pinterpolation inter);

 /* ------------------------------------------------------------
 * Set up polynomial interpolation and evaluate 
 * ------------------------------------------------------------ */

/* Set up chebyshev interpolation points on [a,b]

   'inter' is used to store the interpolation points
   'a' start point of the inerpolation intervall 
   'b' end point of the interpolation intervall */
void
setup_chebyshev_interpolationpoints(pinterpolation inter, double a, double b);

/* Set up aequidistant interpolation points on [a,b] */

void 
setup_aequidistant_interpolationpoints(pinterpolation inter, double a, double b);

/*  Evaluate function f in interpolation points */

void
eval_interpolated_values(pinterpolation inter, function f, void *data);

/* Evaluate Newton divided differences */

void
newton_divided_differences(pinterpolation inter);

/* Use Horner schema to evaluate interpolation polynomial in x */

double
eval_interpolation_polynomial(const pinterpolation inter, double x);


  #endif
