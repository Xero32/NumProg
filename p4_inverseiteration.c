/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 4 - Eigenwertbestimmung (Inverse Iteration) 	 */
/* ------------------------------------------------------------- */
/*	Autoren:	Ove Hansen und Tim Drevelow			 		 */
/*	Versionsnummer:		1.0				 */
/*---------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "miniblas.h"
#include <time.h>
#define EPSILON 1.0e-006


/* ------------------------------------------------------------
 * LR decomposition and solve
 * ------------------------------------------------------------ */
 
static void
tridiag_lrdecomp(ptridiag a){
	//Aufwand 3n
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
	//Aufwand 5n +3n=8n 
	ptridiag h=clone_tridiag(a);
	
	int i;
	int n=h->rows;
	double* l=h->l;
	double* d=h->d;
	double* u=h->u;
	double* x=b->x;
	
	tridiag_lrdecomp(h);
	
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
	
	del_tridiag(h);
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
	
	while(--steps>0){
		double norm=nrm2(n, x->x,1);
		tridiag_lrsolve(a,x);
		scal(n,1.0/norm,x->x,1);
	}
	
	double norm=nrm2(n, x->x,1);
	tridiag_lrsolve(a,x);
	
	pvector y=new_zero_vector(n);
	for(int i=0; i<n; i++){
		y->x[i]=x->x[i];
	}
	
	scal(n,1.0/norm,y->x,1);
	*eigenvalue=dot(n, x->x, 1, y->x, 1)/dot(n, x->x, 1, x->x, 1);
	printf("Steps:%d\n",steps);
	*res = fabs(nrm2(n, y->x, 1)/ *eigenvalue);
}
 
 /* Note: 
   'shift' describs the parameter mu, used to compute (A - mu * I).
   */
 
static void
inverse_iteration_withshift(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){

int n = x->rows,i;


pvector y = new_zero_vector(n);
mvm_tridiag(0, 1, a, x, y);
*eigenvalue = dot(n, x->x, 1, y->x, 1)/dot(n, x->x, 1, x->x, 1);
scal(n,1.0/ fabs(*eigenvalue),y->x,1); //
ptridiag b = clone_tridiag(a);
double* Bd = b->d;
for(i=0 ; i < n; i++){
	Bd[i] -= shift;
}

while(norm2_diff_vector(y,x)> EPSILON && steps-->0){
	double normx = nrm2(n,x->x,1);
	tridiag_lrsolve( b, x);	
	scal(n,1.0/normx,x->x,1);
	clear_vector(y);
	mvm_tridiag(0, 1, a, x, y);
	assert(dot(n, x->x, 1, x->x, 1));
	*eigenvalue = dot(n, x->x, 1, y->x, 1)/dot(n, x->x, 1, x->x, 1);
	scal(n,1.0/fabs(*eigenvalue),y->x,1); //
	//axpy( n, -*eigenvalue, x->x, 1, y->x, 1);
}
	printf("Steps:%d\n",steps);
	*res = norm2_diff_vector(y,x); //
}

 /* Note: 
   In this case 'shift' is only the start parameter mu.
   */

static void
rayleigh_iteration(ptridiag a, pvector x, double shift, int steps, double *eigenvalue, double *res){

int n = x->rows,i;

pvector y = new_zero_vector(n);
mvm_tridiag(0, 1, a, x, y);
*eigenvalue = dot(n, x->x, 1, y->x, 1)/dot(n, x->x, 1, x->x, 1);
scal(n,1.0/ fabs(*eigenvalue),y->x,1); //
/*ptridiag b = clone_tridiag(a);
double* Bd = b->d;
for(i=0 ; i < n; i++){
	Bd[i] -= shift;
}*/

while(norm2_diff_vector(y,x)> EPSILON && steps-->0){
	ptridiag b = clone_tridiag(a);
	double* Bd = b->d;
	for(i=0 ; i < n; i++){
		Bd[i] -= shift;
	}
	double normx = nrm2(n,x->x,1);
	tridiag_lrsolve( b, x);
	del_tridiag(b);
	scal(n,1.0/normx,x->x,1);
	clear_vector(y);
	mvm_tridiag(0, 1, a, x, y);
	*eigenvalue = dot(n, x->x, 1, y->x, 1)/dot(n, x->x, 1, x->x, 1);
	scal(n,1.0/fabs(*eigenvalue),y->x,1);
	//ptridiag b = clone_tridiag(a);
	shift = *eigenvalue;
	//double* Bd = b->d;
	//for(i=0 ; i < n; i++){
	//Bd[i] -= shift;
	//}
}
	printf("Steps:%d\n",steps);
	*res = norm2_diff_vector(y,x);
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
  inverse_iteration(a, x, 20, &lambda, &norm);

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
  mu = 999.0;
  inverse_iteration_withshift(a, x, mu, 20, &lambda, &norm);

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
  mu = 999.0;
  rayleigh_iteration(a, x, mu, 1000, &lambda, &norm);

  printf("  Eigenvalue %g\n"
	 "  Residual norm %e\n",
	 lambda, norm);

  del_vector(x);
  del_tridiag(a);	 


  return EXIT_SUCCESS;
}
