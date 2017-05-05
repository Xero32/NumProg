/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 3 - Block-LR-Zerlegung und Vergleich	 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 				 		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* Additional libraries */
#include "basic.h"		/* basic types and time measurement */
#include "miniblas.h"		/* our simple BLAS version */
#include "matrix.h"		/* matrix functions */


/* older versions e.g. used in comparsion*/
/* ************************************************* */

/* P1 */
static void
lr_decomp(pmatrix a){

    int i,j,k;
    int lda = a->ld;
    int n = a->rows;
    double *aa = a->a;

    for (k = 0; k < n; k++){
        for(i = k+1; i < n; i++){
            aa[i+k*lda] /= aa[k+k*lda];  //assert; nicht durch 0 teilen
        }
        for(i = k+1; i < n; i++){
            for(j = k+1; j < n; j++){
                aa[i+j*lda] -= aa[i + k*lda] * aa[k + j*lda];
            }
        }
    }
}

/* P2 (with slightly changed name) */

static void
lr_decomp_blas(pmatrix a){

    int lda = a->ld;
    double *aa = a->a;
    int k;
    for(k = 0; k < lda; k++){
        assert(aa[k+k*lda]);
        scal(lda-k-1, 1.0/aa[k+k*lda], a->a+k+1+k*lda, 1);
        ger(lda-k-1, lda-k-1, -1.0,
            a->a+(k+1)+k*lda, 1, a->a+k+(k+1)*lda, lda,
            a->a+(k+1)+(k+1)*lda, lda);
    }
}

static void //TODO
block_lsolve(int n, int m, const real *L, int ldL, real *B, int ldB){
for(k=0; k<n; k++)
    int k;
    ger(n-k-1, m, -1.0, L+(k+1)+k*ldL, 1, B+k, ldB, B+(k+1), ldB);
}

static void
block_rsolve_trans(int n, int m, const real *R, int ldR, real *B, int ldB){
    int k;
    for(k=0; k<n; k++) {
        scal(m, 1.0/R[k+k*ldR], B+k*ldB, 1);
        ger(m, n-k-1, -1.0, B+k*ldB, 1, R+k+(k+1)*ldR, ldR, B+(k+1)*ldB, ldB);
    }
}

static void
blocklr_decomp(pmatrix a, int m){
    int i, j, k;
    int lda = a->ld;
    int n = a-> rows;
    double *A = a->a;
    int oi, oj, ok, ni, nj, nk;
    for(k=0; k<m; k++) {
        ok = n * k / m; nk = n * (k+1) / m - ok;
        init_sub_matrix(pmatrix asub, a, nk, ok, nk, ok); //(asub, a, rows, roff, cols, coff)
        lr_decomp_blas(asub);     //(nk, A+ok+ok*ldA, ldA);
        for(j=k+1; j<m; j++) {
            oj = n * j / m; nj = n * (j+1) / m - oj;
            block_lsolve(nk, nj, A+ok+ok*ldA, ldA, A+ok+oj*ldA, ldA);
            block_rsolve_trans(nk, nj, A+ok+ok*ldA, ldA, A+oj+ok*ldA, ldA);
        }
        for(j=k+1; j<m; j++) {
        oj = n * j / m; nj = n * (j+1) / m - oj;
            for(i=k+1; i<m; i++) {
            oi = n * i / m; ni = n * (i+1) / m - oi;
            gemm(false, false, ni, nj, nk, -1.0,
            A+oi+ok*ldA, ldA, A+ok+oj*ldA, ldA, 1.0, A+oi+oj*ldA, ldA);
            }
        }
    }
}


/* ============================================================
 * Main program
 * ============================================================ */

int 
main(void){

  int n;
  pstopwatch sw;
  pmatrix A;
  real time;
  int m;

  n = 2000;					/* matrix dimension */
  m = 100;					/* number of matrix parts */


  
  /* ------------------------------------------------------------
   * Block-LR decomposition
   * ------------------------------------------------------------ */
   init_sub_matrix 
  
  
  
  /* ------------------------------------------------------------
   * 'only' BLAS-LR decomposition
   * ------------------------------------------------------------ */
  /* ------------------------------------------------------------
   * first version of LR decomposition
   * ------------------------------------------------------------ */
  
 
  /* cleaning up */
  del_stopwatch(sw);

  return EXIT_SUCCESS;
}
