/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 3 - Block-LR-Zerlegung und Vergleich	 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marvin Becker, Marko Hollm			 		 */
/*	Versionsnummer:	1					 */
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
            assert(aa[k+k*lda]);
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
    int k;
    for(k=0; k<n; k++){
        ger(n-k-1, m, -1.0, L+(k+1)+k*ldL, 1, B+k, ldB, B+(k+1), ldB);
    }
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
    int ldA = a->ld;
    int n = a-> rows;
    double *A = a->a;
    int oi, oj, ok, ni, nj, nk;
    
    for(k=0; k<m; k++) {
        ok = n * k / m; nk = n * (k+1) / m - ok;
        pmatrix asub = new_matrix(nk,nk);
        asub = init_sub_matrix(asub, a, nk, ok, nk, ok);
        
        lr_decomp_blas(asub);    
        double *AS = asub->a;

        for(j = ok; j < ok+nk; j++){
            for(i = ok; i < ok+nk; i++){
                int u = i-ok;
                int v = j-ok;
                A[i+j*ldA] = AS[u+v*nk];
            }
        }      
        
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
//      del_matrix(asub);        
    }
}


/* ============================================================
 * Main program
 * ============================================================ */

int 
main(void){

  int n;
  pmatrix A,B,C;
  real time1,time2,time3;
  int m;
  int max;
  FILE *f = NULL;
  n = 2000;     /* matrix dimension */
  m = 100;
  int ctr = 1;
  max = 0;
				 
// switch(n){
//     case 1000: f = fopen("data1.dat","w"); max = 10; break;
//     case 2000: f = fopen("data2.dat","w"); max = 11; break;
//     case 3000: f = fopen("data3.dat","w"); max = 14; break;
//     case 4000: f = fopen("data4.dat","w"); max = 17; break;
//     case 5000: f = fopen("data5.dat","w"); max = 18; break;
//     case 6000: f = fopen("data6.dat","w"); max = 18; break;
// }
// 
//     REPEAT:
//   switch(ctr){   /* number of matrix parts */
// //       case 1: m = 16; break;
//       case 2: m = n/100; break;
//       case 3: m = 32; break;
//       case 4: m = 50; break;
//       case 5: m = 64; break;
//       case 6: m = 100; break;
//       case 7: m = 128; break;
//       case 8: m = 150; break;
//       case 9: m = 200; break;
//       case 10: m = 250; break;
//       case 11: m = 300; break;
//       case 12: m = 350; break;
//       case 13: m = 400; break;
//       case 14: m = 450; break;
//       case 15: m = 500; break;
//       case 16: m = 550; break;
//       case 17: m = 600; break;
//       case 18: m = 650; break;
//   }
  pstopwatch sw = new_stopwatch();
  /* ------------------------------------------------------------
   * Block-LR decomposition
   * ------------------------------------------------------------ */
  C = new_diaghilbert_matrix(n);
  B = new_diaghilbert_matrix(n);
  A = new_diaghilbert_matrix(n);

  start_stopwatch(sw);
  blocklr_decomp(A,m);  
//   printf("Block decomp:\n");
//   print_matrix(A);
  printf("Duration of Block decomp: %f\n",time1 = stop_stopwatch(sw));
  
  /* ------------------------------------------------------------
   * 'only' BLAS-LR decomposition
   * ------------------------------------------------------------ */
    start_stopwatch(sw);
    lr_decomp_blas(B);
//     printf("BLAS decomp:\n");
//     print_matrix(B);
      
    printf("Duration of BLAS decomp: %f\n",time2 = stop_stopwatch(sw));
  
  /* ------------------------------------------------------------
   * first version of LR decomposition
   * ------------------------------------------------------------ */
//   if(n < 4001){
//       if(m == 16){
    start_stopwatch(sw);
    lr_decomp(C);
//     printf("Basic decomp:\n");
//     print_matrix(C);
    printf("Duration of basic decomp: %f\n",time3 = stop_stopwatch(sw));
//     } 
//   }
  /* ------------------------------------------------------------
   * test functioning
   * ------------------------------------------------------------ */    
    double *aa = A->a;
    double *ba = B->a;
    double err1 = 0.0;
    for(int j = 0; j < n; j++){
        for(int i = 0; i < n; i++){
            err1 = (err1 < fabs(aa[i+j*n] - ba[i+j*n])) ? fabs(aa[i+j*n] - ba[i+j*n]) : err1;
//             printf("index: %d\t errval: %f\t, relative err: %f\n",i+j*n,fabs(aa[i+j*n] - ba[i+j*n]),fabs(aa[i+j*n] - ba[i+j*n])/ba[i+j*n]);
        }
    }
    printf("Error Value: %f\n",err1);
    
    
  /* cleaning up */
  del_stopwatch(sw);
  del_matrix(A);
  del_matrix(B);
  del_matrix(C);

  if(f){
      fprintf(f,"%d\t%f\t%f\t%f\n",m,time1,time2,time3);
  }

  if(ctr <= max){
      ctr++;
//       goto REPEAT;
  }
  
  if(f)  fclose(f);
  return EXIT_SUCCESS;
}
