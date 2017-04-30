/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 2 - LR-Zerlegung and inversion inkl. BLAS	 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker         	 */
/*	Versionsnummer:	0.5					 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* Additional libraries */
#include "basic.h"		/* basic types and time measurement */
#include "miniblas.h"		/* our simple BLAS Version */
#include "matrix.h"		/* matrix functions */


/* ------------------------------------------------------------
 * LR decomposition using BLAS
 * ------------------------------------------------------------ */

static void //TODO
lrdecomp(pmatrix a){

    int lda = a->ld;
    double *aa = a->a;
    int k;
    for(k = 0; k < lda; k++){
        assert(aa[k+k*lda]);
        scal(lda-k-1, 1.0/aa[k+k*lda], aa+k+1+k*lda, 1);
        ger(lda-k-1, lda-k-1, -1.0,
            aa+(k+1)+k*lda, 1, aa+k+(k+1)*lda, lda,
            aa+(k+1)+(k+1)*lda, lda);
    }
}


static void //TODO
lowersolve_matrix(int unit, const pmatrix a, pvector b){
    int k;
    int lda = a->ld;
    double *aa = a->a;
    double *bb = b->x;
    
    for(k=0; k<unit; k++){
        axpy(unit-k-1, -bb[k], aa+(k+1)+k*lda, 1, bb+k+1, 1);
    }
}


static void //TODO
uppersolve_matrix(int unit, const pmatrix a, pvector b){
    int k;
    int lda = a->ld;
    double *aa = a->a;
    double *bb = b->x;
    for(k=unit; k-->0;){
        assert(aa[k+k*lda]);
        bb[k] /= aa[k+k*lda];
        axpy(k, -bb[k], aa+k*lda, 1, bb, 1);
    }
    
}


static pmatrix //TODO
invert(const pmatrix a){
    int n = a->rows;
    int lda = a->ld;
    pmatrix z = new_zero_matrix(n,n);
    
    double *aa = a->a;
    double *zz = z->a;
    
    
    for(int j = n-1; j--> 0 ;){
        for(int i = n; i--> 0;){
            zz[i+j*lda] = aa[i+j*lda];  //z = lower matrix
            aa[i+j*lda] = 0.0;
        }
    }
    for(int k = 0; k < n; k++){
        zz[k+k*lda] = 1.0;
    }
    
    for(int i = 0; i < n; i++){
        pvector e = new_zero_vector(n);
        double *ee = e->x;
        ee[i] = 1.0;
        uppersolve_matrix(n,a,e);
        for(int j = 0; j < n; j++){
            aa[j+i*lda] = ee[j];
        }   
        del_vector(e);
    }
    
    for(int i = 0; i < n; i++){
        pvector f = new_zero_vector(n);
        double *ff = f->x;
        ff[i] = 1.0;
        lowersolve_matrix(n,z,f);
        for(int j = 0; j < n; j++){
            zz[j+i*lda] = ff[j];
        }
        del_vector(f);
    }
    for (int j = n-1; j--> 0;){
        for(int i = n; i--> 0;){
            aa[i+j*lda] = zz[i+j*lda];
        }
    }
        
    return z;
}

/* ============================================================
 * Main program
 * ============================================================ */

int
main(void){

  pmatrix a, inva;
  pvector b, x, xrev;
  double err, err_abs, norm_x;
  int rows;
  int i;
  
  rows = 12;
  
  /* ------------------------------------------------------------
   * Hilbert matrix, LR decomposition
   * ------------------------------------------------------------ */
  b  = new_zero_vector(rows);
  x  = new_zero_vector(rows);
  xrev  = new_vector(rows);

  printf("Testing Hilbert matrix, LR decomposition\n");

  a  = new_hilbert_matrix(rows);

  for(i = 0; i < rows; i++){
    xrev->x[i] = 1.0 / (1.0 + i);
  }

  gemv(false, a->rows, a->cols, 1.0, a->a, a->ld, xrev->x, 1, b->x, 1);
  printf("matrix a\n");
  print_matrix(a); //TODO
  lrdecomp(a);
  printf("decomp matrix a\n");
  print_matrix(a); //TODO
  inva = invert(a);
  gemv(false, a->rows, a->cols, 1.0, inva->a, a->ld, b->x, 1, x->x, 1);
  printf("inverted matrix\n");
  print_matrix(inva); //TODO
  err_abs = 0.0;
  norm_x  = 0.0;
  for(i = 0; i < rows; i++) {
    err = fabs(xrev->x[i] - x->x[i]);
    err_abs = (err_abs < err ? err : err_abs);
    norm_x = (norm_x < fabs(xrev->x[i]) ? fabs(xrev->x[i]) : norm_x);
  }
  printf("  Absolute max error: %.3e\n", err_abs);
  printf("  Relative max error: %.3e\n", err_abs/norm_x);

  del_matrix(a);
  //del_matrix(inva); //TODO
  del_vector(x);
  del_vector(xrev);
  del_vector(b);

  /* ------------------------------------------------------------
   * Diagonally dominant matrix, QR decomposition
   * ------------------------------------------------------------ */

  b  = new_zero_vector(rows);
  x  = new_zero_vector(rows);
  xrev  = new_vector(rows);

  printf("Testing diagonally dominant Hilbert matrix, LR decomposition\n");

  a  = new_diaghilbert_matrix(rows);

  for(i = 0; i < rows; i++){
    xrev->x[i] = 1.0 / (1.0 + i);
  }

  gemv(false, a->rows, a->cols, 1.0, a->a, a->ld, xrev->x, 1, b->x, 1);
  printf("matrix a\n");
  print_matrix(a);
  lrdecomp(a);
  printf("decomp matrix a\n");
  print_matrix(a);
  inva = invert(a);
  gemv(false, a->rows, a->cols, 1.0, inva->a, a->ld, b->x, 1, x->x, 1);
  
  err_abs = 0.0;
  norm_x  = 0.0;
  for(i = 0; i < rows; i++) {
    err = fabs(xrev->x[i] - x->x[i]);
    err_abs = (err_abs < err ? err : err_abs);
    norm_x = (norm_x < fabs(xrev->x[i]) ? fabs(xrev->x[i]) : norm_x);
  }
  printf("  Absolute max error: %.3e\n", err_abs);
  printf("  Relative max error: %.3e\n", err_abs/norm_x);
    
  //lowersolve_matrix(rows,a,b);
  //uppersolve_matrix(rows,a,b);
  del_matrix(a);
  //del_matrix(inva); //TODO
  del_vector(x);
  del_vector(xrev);
  del_vector(b);

  
  return EXIT_SUCCESS;
}
