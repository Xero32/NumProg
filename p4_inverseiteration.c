/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 4 - Eigenwertbestimmung (Inverse Iteration) 	 */
/* ------------------------------------------------------------- */
/*	Autoren:				 	Ove Hansen, Tim Drevelow	 */
/*	Versionsnummer:				2.0		 */
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
	int i;
	int n=a->rows;
	double* l=a->l;
	double* d=a->d;
	double* u=a->u;
	for(i=0; i<n-1; i++){
		assert(d[i]);
		l[i]/=d[i];
		d[i+1]-=l[i]*u[i];
	}
  
}

static void
tridiag_lrsolve(const ptridiag a, pvector b){
	
	//create clone of a, because should not be changed
	ptridiag h=clone_tridiag(a);
	assert(a->rows==b->rows);
	
	int i;
	int n=h->rows;
	double* l=a->l;
	double* d=a->d;
	double* u=a->u;
	double* x=b->x;
	
	
	//solve L*y=b
	for(i=1; i<n; i++){
		x[i]-=l[i-1]*x[i-1];
	}
	
	//solve R*x=y
	assert(d[n-1]);
	x[n-1]/=d[n-1];
	for(i=n-1; i>0; i--){
		assert(d[i-1]);
		x[i-1]=(x[i-1]-x[i]*u[i-1])/d[i-1];
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

static void
inverse_iteration(ptridiag a, pvector x, int steps, double *eigenvalue, double *res){
	
	int n=x->rows;
	double nrmx;
	ptridiag b;
	
	b=clone_tridiag(a);
	tridiag_lrdecomp(b);
	
	while(steps-->0){
		nrmx = nrm2(n, x->x,1);
		assert(nrmx);
		scal(n,1.0/nrmx,x->x,1);
		tridiag_lrsolve(b,x);
	}
	
	nrmx = nrm2(n, x->x,1);
	assert(nrmx);
	scal(n,1.0/nrmx,x->x,1);
	pvector y=new_zero_vector(n);
	mvm_tridiag(0, 1, a, x, y);
	assert(dot(n, x->x, 1, x->x, 1));
	*eigenvalue=dot(n, x->x, 1, y->x, 1);
	axpy(n, -*eigenvalue,x->x, 1,y->x, 1);
	*res=nrm2(n,y->x,1);
	
	
}
 
 /* Note: 
   'shift' describs the parameter mu, used to compute (A - mu * I).
   */
static void
inverse_iteration_withshift(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){


	int n=x->rows;
	double nrmx;
	int i;
	ptridiag c=clone_tridiag(a);
	for(i=0; i<n; i++) c->d[i]-=shift;
	tridiag_lrdecomp(c);
	
	while(steps-->0){
		nrmx = nrm2(n, x->x,1);
		assert(nrmx);
		scal(n,1.0/nrmx,x->x,1);
		tridiag_lrsolve(c,x);
	}
	
	nrmx = nrm2(n, x->x,1);
	assert(nrmx);
	scal(n,1.0/nrmx,x->x,1);
	pvector y=new_zero_vector(n);
	mvm_tridiag(0, 1, a, x, y);
	assert(dot(n, x->x, 1, x->x, 1));
	*eigenvalue=dot(n, x->x, 1, y->x, 1);
	axpy(n, -*eigenvalue,x->x, 1,y->x, 1);
	*res=nrm2(n,y->x,1);
	
	del_tridiag(c);

}

 /* Note: 
   In this case 'shift' is only the start parameter mu.
   */

static void
rayleigh_iteration(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){

	int i;
	int n=x->rows;
	double nrmx;
	pvector y=new_vector(n);
	ptridiag c=clone_tridiag(a);
	
	for(i=0; i<n; i++) c->d[i]-=shift;
	nrmx = nrm2(n, x->x,1);
	assert(nrmx);
	scal(n,1.0/nrmx,x->x,1);
	tridiag_lrdecomp(c);
	tridiag_lrsolve(c,x);
	
	
	while(steps-->0){
		
		nrmx = nrm2(n, x->x,1);
		assert(nrmx);
		scal(n,1.0/nrmx,x->x,1);		
		clear_vector(y);
		mvm_tridiag(0, 1, a, x, y);
		assert(dot(n, x->x, 1, x->x, 1));
		shift=dot(n, x->x, 1, y->x, 1);
		
		copy_tridiag(a,c);
	    for(i=0; i<n; i++) c->d[i]-=shift;
		
		tridiag_lrdecomp(c);
		tridiag_lrsolve(c,x);
		
	}
	del_tridiag(c);
	
	nrmx = nrm2(n, x->x,1);
	assert(nrmx);
	scal(n,1.0/nrmx,x->x,1);	
	clear_vector(y);
	mvm_tridiag(0, 1, a, x, y);
	assert(dot(n, x->x, 1, x->x, 1));
	*eigenvalue=dot(n, x->x, 1, y->x, 1);
	axpy(n, -*eigenvalue,x->x, 1,y->x, 1);
	*res=nrm2(n,y->x,1);
	
	
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

  n = 1000;

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
  inverse_iteration(a, x, 100, &lambda, &norm);

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
  inverse_iteration_withshift(a, x, mu, 100, &lambda, &norm);

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
  rayleigh_iteration(a, x, mu, 100, &lambda, &norm);

  printf("  Eigenvalue %g\n"
	 "  Residual norm %e\n",
	 lambda, norm);

  del_vector(x);
  del_tridiag(a);	 


  return EXIT_SUCCESS;
}
