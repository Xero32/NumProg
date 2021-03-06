/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 4 - Eigenwertbestimmung (Inverse Iteration) 	 */
/* ------------------------------------------------------------- */
/*	Autoren:	Marko Hollm, Marvin Becker		 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "miniblas.h"
#include <time.h>

#define EPSILON 0.0000000000001
#define PRNTTD(Lm,Dm,Um,i) printf("L: %f\tD: %f\tU: %f\n",Lm[i],Dm[i],Um[i]);
#define PRNTVC(bb,j) printf("%f, ",bb[j]);
#define FOR for(i=0; i < n; i++)
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
        D[k+1] -= L[k]*U[k];                        // 3n steps
    }
}

static void
tridiag_lrsolve(const ptridiag a, pvector b){
    assert(a->rows == b->rows);
    ptridiag f = clone_tridiag(a);
    double* L = f->l;
    double* D = f->d;
    double* U = f->u;
    double* B = b->x;
    int n = a->rows;
    tridiag_lrdecomp(f);
    /* L_solve, Ly=b */
    for(int k = 0; k < n-1; k++){
        assert(L[k] && B[k]);
        B[k+1] -= L[k]*B[k];                        // 2n steps
    }
    /* R_solve, Rx=y */
    int m = n-1;
    assert(D[m]);
    B[m] /= D[m];                       //for(k = n; k-- > 0;) loop from n-1 to 0
    for(int k = m; k-- > 0;){
        assert(D[k]);
        B[k] = (B[k] - B[k+1] * U[k]) / D[k];       // 3n steps
    }
    del_tridiag(f);
}
                                                            // => 8n steps in fct lrsolve
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
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax

   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x
   
   assert(*eigenvalue);
   scal(n, 1/fabs(*eigenvalue), y->x, 1);
   while(norm2_diff_vector(x,y) > EPSILON && steps-- > 0){                      // while ||y - lam x|| > eps lam   
        normx = nrm2(n,x->x,1); 
        assert(normx);
        scal(n, 1.0/normx, x->x, 1);                                            // x <- x' / ||x||; first part done in line above; only scaling with inverse of normx done here
        tridiag_lrsolve(a,x);                                                   // A x = x';
        y = new_zero_vector(n);                                                 // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1)); // would assert(x->x); be enough? //TODO
        *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x

        assert(*eigenvalue);
        scal(n, 1.0/fabs(*eigenvalue), y->x, 1);                                // scaling for termination condition
    }   
   *res = (norm2_diff_vector(x,y)) * *eigenvalue;                               // rescaling to get residual norm
   del_vector(y);  
   printf("steps: %d\n",steps);
}

 /* Note: 
   'shift' describes the parameter mu, used to compute (A - mu * I).
   */
 
static void
inverse_iteration_withshift(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){      
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax
   
   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x
   
   ptridiag b = clone_tridiag(a);
   double *Bd = b->d;
   for(int i = 0; i< n; i++){
        Bd[i] -= shift;                                                         // B <- A - my I
   }
   
   assert(*eigenvalue);
   scal(n, 1.0/fabs(*eigenvalue), y->x, 1);

   while(norm2_diff_vector(x,y) > EPSILON && steps-- > 0){                      // while ||y/lam - x|| > eps

        normx = nrm2(n,x->x,1);                                                 
        assert(normx);
        tridiag_lrsolve(b,x); //lrdecomp already happens in lrsolve             // B x' = (A - my I) x' = x; write: x <- x'
        scal(n, 1.0/normx, x->x, 1);                                            // x <- x' / ||x||; first part done above; only scaling with inverse of normx done here
        
        clear_vector(y);                                                        // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1));
        *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x

        assert(*eigenvalue);
        scal(n, 1.0/fabs(*eigenvalue), y->x, 1);
        
   }   
   *res = (norm2_diff_vector(x,y)) * fabs(*eigenvalue);
   del_vector(y);
   printf("steps: %d\n",steps);
}

 /* Note: 
   In this case 'shift' is only the start parameter mu.
   */

static void
rayleigh_iteration(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){        
   assert(a->rows == x->rows);  
   int n = x->rows;
   double normx;
   pvector y = new_zero_vector(n);
   
   mvm_tridiag(0, 1, a, x, y);                                                  // y <- Ax
   
   assert(dot(n, x->x, 1, x->x, 1));
   *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);           // lam <- x*y/x*x
   assert(*eigenvalue);
   scal(n, 1.0/fabs(*eigenvalue), y->x, 1);

   while(norm2_diff_vector(x,y) > EPSILON && steps-- > 0){                      // while ||y/lam - x|| > eps
        ptridiag b = clone_tridiag(a);
        double *Bd = b->d;
        for(int i = 0; i< n; i++){
            Bd[i] -= shift;                                                     // B <- A - my I
        }
        tridiag_lrsolve(b,x); //lrdecomp already happens in lrsolve             // B x' = (A - my I) x' = x; write: x <- x'
        
        normx = nrm2(n,x->x,1);                                                 
        assert(normx);
        scal(n, 1.0/normx, x->x, 1);                                              // x <- x' / ||x||; first part done above; only scaling with inverse of normx done here
        
        clear_vector(y);                                                        // set y as zero-vector again
        mvm_tridiag(0, 1, a, x, y);                                             // y <- Ax; input zero-vec which to write to
        assert(dot(n, x->x, 1, x->x, 1));
        shift = *eigenvalue = dot(n, x->x, 1, y->x, 1) / dot(n, x->x, 1, x->x, 1);      // lam <- x*y/x*x
        assert(*eigenvalue);
        scal(n, 1.0/fabs(*eigenvalue), y->x, 1);
        del_tridiag(b);     //efficient?
   }
   *res = (norm2_diff_vector(x,y)) * fabs(*eigenvalue);
   del_vector(y);
   printf("steps: %d\n",steps);
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

//   double *U,*D,*L;
//   U = a->u;
//   D = a->d;
//   L = a->l;
//   for(int i = 0; i < n; i++){
//         printf("L: %f\tD: %f\tU: %f\n",L[i],D[i],U[i]);
//   }
  
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

  /* ------------------------------------------------------------
   * Testing LR solve
   * ------------------------------------------------------------ */
// ptridiag M = new_threepointstencil(n); 
// pvector b = new_zero_vector(n);
// double *bb = b->x;
// int i;
// FOR{
//     bb[i] = 10*(i+1);
// } 
// double *Um,*Dm,*Lm;
//   Um = M->u;
//   Dm = M->d;
//   Lm = M->l;
// FOR{
//     PRNTTD(Lm,Dm,Um,i);
// }
// 
// FOR{
//     PRNTVC(bb,i);
// }
// printf("\n");
// 
// tridiag_lrsolve(M,b);
// 
// for(int j = 0; j < n; j++){
//     printf("%f, ",bb[j]);
// }
// printf("\n");
//   
  
  return EXIT_SUCCESS;
}

/*  EIGENVALUES calculated with octave for n = 10:
 * 
    9.8027
    38.4166
    83.5237
   141.4696
   207.5598
   276.4402
   342.5304
   400.4763
   445.5834
   474.1973
   */