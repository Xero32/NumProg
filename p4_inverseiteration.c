/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 4 - Eigenwertbestimmung (Inverse Iteration) 	 */
/* ------------------------------------------------------------- */
/*	Autoren:				 		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "miniblas.h"
#include <time.h>


/* ------------------------------------------------------------
 * LR decomposition and solve
 * ------------------------------------------------------------ */

static void
tridiag_lrdecomp(ptridiag a){
    
    double* L = a->l;
    double* D = a->d;
    double* U = a->u;
    int n = a->rows;
    
    for(int k = 0; k < n-1; k++){
        assert(D[k]);
        L[k] /= D[k];
        D[k+1] -= L[k]*U[k];
    }
}

static void
tridiag_lrsolve(const ptridiag a, pvector b){
	
    double* L = a->l;
    double* D = a->d;
    double* U = a->u;
    double* B = b->x;
    int n = a->rows;
    tridiag_lrdecomp(a);
    /* L_solve, Ly=b */
    for(int k = 0; k < n-1; k++){
        assert(L[k] && B[k]);
        B[k+1] /= L[k]*B[k];
    }
    /* R_solve, Rx=y */
    assert(D[n-1]);
    B[n-1] /= D[n-1];                       //for(k = n; k-- > 0;) loop from n-1 to 0
    for(int k = n-1; k-- > 1;){
        assert(D[k-1]);
        B[k-1] = (B[k-1] - B[k] * U[k-1]) / D[k-1];
    }
}

/* ------------------------------------------------------------
 * inverse iteration
 * ------------------------------------------------------------ */

/* Note: 
   'steps' is used to determine the maximum number
   of steps, which should be used to the iteration step.
   'eigenvalue' will be overwritten by the computed 
    eigenvalue.
   'res' is used to determine the error, measured in the 
    euclidean norm.
   */

static void /* copied 'inverse_iteration_withshift' but removed all appearances of shift and auxiliary matrix B */
inverse_iteration(ptridiag a, pvector x, int steps, double *eigenvalue, double *res){
//function may need some further assessment //TODO
   int ctr = 0;                                                                 // counter to keep track of iteration step
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax
   
   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x

   axpy(n, -(*eigenvalue), x->x, 1, y->x, 1);                                   // y <- alpha x + y    with alpha = -lam   =>   y <- y - lam x
   
   while(nrm2(n, y->x, 1) > *res * fabs(*eigenvalue) && ctr < steps){           // while ||y - lam x|| > eps lam   
        /* alternatively one could try something like this: 
         *nrm(y - lam x) > eps lam <=> nrm(y/lam - x) > eps */

        normx = nrm2(n,x->x,1);                                                 
        assert(normx);
        scal(n, 1/normx, x->x, 1);                                              // x <- x' / ||x||; first part done in line above; only scaling with inverse of normx done here

        y = new_zero_vector(n);                                                 // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1)); // would assert(x->x); be enough? //TODO
        *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x
        
        ctr++;
        assert(*eigenvalue);
        *res = nrm2(n, y->x, 1) / *eigenvalue;    //TODO how to use *res?
   }   
   del_vector(y);
}
 
 /* Note: 
   'shift' describes the parameter mu, used to compute (A - mu * I).
   */
 
static void
inverse_iteration_withshift(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){
   int ctr = 0;                                                                 // counter to keep track of iteration step        
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax
   
   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x
   
   axpy(n, -(*eigenvalue), x->x, 1, y->x, 1);                                   // y <- alpha x + y    with alpha = -lam   =>   y <- y - lam x
   ptridiag b = clone_tridiag(a);
   
   double *Bd = b->d;
   for(int i = 0; i< n; i++){
        Bd[i] -= shift;                                                         // B <- A - my I
   }

   while(nrm2(n, y->x, 1) > *res * fabs(*eigenvalue) && ctr < steps){                          // while ||y - lam x|| > eps lam
        /* alternatively one could try something like this: 
         * nrm(y - lam x) > eps lam <=> nrm(y/lam - x) > eps */
        
        normx = nrm2(n,x->x,1);                                                 
        assert(normx);
        tridiag_lrsolve(b,x); //lrdecomp already happens in lrsolve             // B x' = (A - my I) x' = x; write: x' <- x
        scal(n, 1/normx, x->x, 1);                                              // x <- x' / ||x||; first part done in line above; only scaling with inverse of normx done here

        y = new_zero_vector(n);                                                 // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1));
        *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x
        
        ctr++;
        assert(*eigenvalue);
        *res = nrm2(n, y->x, 1) / *eigenvalue;
   }   
   del_vector(y);
}

 /* Note: 
   In this case 'shift' is only the start parameter mu.
   */

static void
rayleigh_iteration(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){
    //TODO
   int ctr = 0;                                                                 // counter to keep track of iteration step
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax
   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x
   
   axpy(n, -(*eigenvalue), x->x, 1, y->x, 1);                                   // y <- alpha x + y    with alpha = -lam   =>   y <- y - lam x
   ptridiag b = clone_tridiag(a);
//    del_vector(y);
   double *Bd = b->d;
   

   while(nrm2(n, x->x, 1) > *res * fabs(*eigenvalue) && ctr < steps){                          // while ||y - lam x|| > eps lam
       /* alternatively one could try something like this: 
        * nrm(y - lam x) > eps lam <=> nrm(y/lam - x) > eps */
       shift = *eigenvalue;                                                     // new in Rayleigh-Iteration: shift varies each time
       for(int i = 0; i< n; i++){                                               // therefore we need to solve " B x' = (A - my I) x' = x " in each iteration    
            Bd[i] -= shift;                                                     // B <- A - my I
        }
               
        normx = nrm2(n,x->x,1);                                                 
        assert(normx);
        tridiag_lrsolve(b,x); //lrdecomp already happens in lrsolve             // B x' = (A - my I) x' = x; write: x' <- x
        scal(n, 1/normx, x->x, 1);                                              // x <- x' / ||x||; first part done in line above; only scaling with inverse of normx done here

        y = new_zero_vector(n);                                                 // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1));
        *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x
        
        ctr++;
        assert(*eigenvalue);
        *res = nrm2(n, y->x, 1) / *eigenvalue;
   } 
   del_vector(y);
}
 
/* ============================================================
 * Main program
 * ============================================================ */

int
main(void){

  tridiag *a;
  vector *x;
  double norm, lambda, mu;
  int n;
  time_t t;

  n = 100;

  time(&t);
  srand((unsigned int) t);

  a = new_threepointstencil(n);
  x = new_vector(n);

  random_vector(x);	

  /* ------------------------------------------------------------
   * Inverse iteration
   * ------------------------------------------------------------ */

  printf("Inverse iteration\n");
  lambda = 0.0;
  norm = 0.0;
  inverse_iteration(a, x, 10, &lambda, &norm);

  printf("  Eigenvalue %g\n"
	 "  Residual norm %e\n",
	 lambda, norm);
	 
 /* ------------------------------------------------------------
   * Inverse iteration with shift
   * ------------------------------------------------------------ */
  random_vector(x); 

  printf("Inverse iteration with shift\n");
  lambda = 0.0;
  norm = 0.0;
  mu = 7.0;
  inverse_iteration_withshift(a, x, mu, 10, &lambda, &norm);

  printf("  Eigenvalue %g\n"
	 "  Residual norm %e\n",
	 lambda, norm);

  /* ------------------------------------------------------------
   * Rayleigh iteration
   * ------------------------------------------------------------ */

  printf("Rayleigh iteration\n");
 
  random_vector(x);
  lambda = 0.0;
  norm = 0.0;
  mu = 7.0;
  rayleigh_iteration(a, x, mu, 10, &lambda, &norm);

  printf("  Eigenvalue %g\n"
	 "  Residual norm %e\n",
	 lambda, norm);

  del_vector(x);
  del_tridiag(a);	 


  return EXIT_SUCCESS;
}
   // nrm(y - lam x) > eps lam <=> nrm(y/lam - x) > eps;