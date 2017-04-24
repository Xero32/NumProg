/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 		Serie 1 - Inversion mit LR-Zerlegung	 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker		 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/* -------------------------------------------------
   Matrix struct and auxiliary functions 
   ------------------------------------------------- */


typedef struct {
  int rows;
  int cols;
  int ld;

  double *a;
} matrix;

typedef matrix *pmatrix;

static pmatrix
new_matrix(int rows, int cols){

  pmatrix a;

  a = (pmatrix) malloc(sizeof(matrix));
  a->a = (double *) malloc(sizeof(double) * rows * cols);
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;
  
  return a;
}

static pmatrix
new_zero_matrix(int rows, int cols){

  pmatrix a;
  int i, j;

  a = (pmatrix) malloc(sizeof(matrix));
  a->a = (double *) malloc(sizeof(double) * rows * cols);
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;
  
  for(i = 0; i < rows; i++){
    for(j = 0; j < cols; j++){
      a->a[i + j *a->ld] = 0.0;
    }
  }

  return a;
}

void
del_matrix(pmatrix a){

  free(a->a);
  free(a);
}


void
print_matrix(pmatrix a){

  int i, j;
  int lda = a->ld;

  printf("Matrix (%d, %d)\n", a->rows, a->cols);

  for(i = 0; i < a->rows; i++){
    printf("( %f ", a->a[i]);
    for(j = 1; j < a->cols; j++){
       printf(", %f ", a->a[i + j*lda]);
    }
    printf(")\n");
  }

}

/* -------------------------------------------------
   Some kind of invertible matrices
   ------------------------------------------------- */

/* Simple 2 times 2 version */
pmatrix
new_2x2_matrix(void){
	
  pmatrix a;
  double *aa;
	
  a = new_matrix(2,2);
  aa = a->a;
	
  aa[0] = 2.0;
  aa[1] = 4.0;
  aa[2] = -1.0;
  aa[3] = 0.0;
	
  return a;	
}

/* Simple 3 times 3 version */
pmatrix
new_3x3_matrix(void){
	
  pmatrix a;
  double *aa;
	
  a = new_matrix(3,3);
  aa = a->a;
	
  aa[0] = 4.0;
  aa[1] = 5.0;
  aa[2] = 3.0;
  aa[3] = 3.0;
  aa[4] = 6.0;
  aa[5] = 2.0;
  aa[6] = 1.0;
  aa[7] = 0.0;
  aa[8] = 2.0;
	
  return a;	
}


/* Simple 4 times 4 version */
pmatrix
new_4x4_matrix(void){

  pmatrix a;
  double *aa;
    
  a = new_matrix(4, 4);
  aa = a->a;  
	
  aa[0] =2.0;
  aa[1] =4.0;
  aa[2] =6.0;
  aa[3] =-2.0;
  aa[4] =-1.0;
  aa[5] =0.0;
  aa[6] =1.0;
  aa[7] =-5.0;
  aa[8] =-3.0;
  aa[9] =-3.0;
  aa[10] =-1.0;
  aa[11] =4.0;
  aa[12] =3.0;
  aa[13] =1.0;
  aa[14] =6.0;
  aa[15] =1.0;	
	
  return a;	
}

/* One ill-conditioned matrix */
pmatrix
new_hilbert_matrix(int rows){

  pmatrix a;
  double *aa;
  int lda;
  int i, j;

  a = new_matrix(rows, rows);
  aa = a->a;
  lda = a->ld;

  for(j = 0; j < rows; j++){
    for(i = 0; i < rows; i++){
      aa[i+j*lda] = 1.0 / (1.0 + i + j);
    }
  }

  return a;
}


/* Finite difference discretization matrix 1D*/
pmatrix
new_3pointstencil(int n){
	
  pmatrix a;
  double *aa;	
  int i;
  int lda;
  
  if(n >1){
  a = new_matrix(n, n);
  aa = a->a;
  lda = a->ld;
	
  aa[0] = 2.0;
  for(i = 1; i < n; i++){
	aa[i + i*lda] = 2.0;
	aa[(i-1) + i*lda] = -1.0;
	aa[i + (i-1)*lda] = -1.0;
  }
  aa[(n-1) + (n-2)*lda] = -1.0;
  }
  else{
	a = new_matrix(n, n);
    a->a[0] = 2;	
  }
  return a;	
}

/* -------------------------------------------------
   Frobenius norm  
   ------------------------------------------------- */


double
normfrob_diff_matrix(pmatrix a){
  
  int i,j;
  double *aa = a->a; 
  double norm = 0.0;  
  int lda = a->ld;
  int n = a->rows;
  double delta;
		
  for(i = 0; i < n; i++){
	for(j = 0; j < n; j++){
		delta = ( i == j ? 1.0 : 0.0);
		norm += fabs((aa[i + j*lda] - delta));
	}  
  }	
	
  return norm;	
}


/* -------------------------------------------------
   MVM, LR-decomposition and so on...  
   ------------------------------------------------- */

/* Very simple matrix multiplication without scalars and transposition */
void
mm(pmatrix a, pmatrix b, pmatrix t){
	
  int i, j, k;
  int lda = a->ld;
  int n = a->rows;
  double *aa = a->a;
  double *ba = b->a;
  double *ta = t->a;  
	
  assert(t->cols == n);
  assert(a->rows == n);  
  assert(a->cols == n);
  assert(b->rows == n);
  assert(b->cols == n);
  
  for(i = 0; i < n; i++){
	for(j = 0; j < n; j++){
		ta[i + j*lda] = 0.0;
	  for(k = 0; k < n; k++){
		ta[i + j*lda] += aa[i + k*lda]* ba[k + j*lda];
	  }
	}
  }  
	
	
}

/* Inplace LR-decomposition without pivot search */
void
lr_decomp(pmatrix a){
    
    int i,j,k;
    int lda = a->ld;
    int n = a->rows;
    double *aa = a->a;

    for (k = 0; k < n; k++){
        for(i = k+1; i < n; i++){
            aa[i+k*lda] /= aa[k+k*lda];
        }
        for(i = k+1; i < n; i++){
            for(j = k+1; j < n; j++){
                aa[i+j*lda] -= aa[i + k*lda] * aa[k + j*lda];
            }
        }
    }
}

/* Inplace inversion of L and R */
void
lr_invert(pmatrix a){

    int i,j,k;
    int lda = a->ld;
    int n = a->rows;
    double *aa = a->a;
    double sum;
    
    //R invert
    for(i = n; i-- > 0;){
        aa[i + i*lda] = 1 / aa[i + i*lda];
        for(j = n; j-- > i+1;){
            sum = 0.0;
            for(k = j; k >= i+1; k--){                   
                sum -= aa[i + k*lda] * aa[k + j*lda];
            }
            aa[i + j*lda] = sum * aa[i + i*lda];
          }
    }
    
    //L invert
    for(i = 1; i < n; i++){
        for(j = 0; j < i; j++){
            aa[i + j*lda] = -aa[i + j*lda];
            for(k = j+1; k < i; k++){
                aa[i + j*lda] -= aa[i + k*lda] * aa[k + j*lda];   
            }
        }
    }
}

/* Inplace multiplication of R^{-1} with L^{-1} */ 
void
lr_mm(pmatrix a){
	
	
  int i, j, k;
  double sum;
  int lda = a->ld;
  int n = a->rows;
  double *aa = a->a;

  assert(a->cols == n);

  for(i = 0; i < n; i++){
	  for(j = 0; j < n; j++){
		if(i <= j){
			sum = aa[i + j*lda];			/* multiplied by l_{jj} = 1*/
		    for(k = j +1; k < n; k++){
				sum += aa[i + k*lda] * aa[k + j*lda];
			}
		}
		else{
			sum = 0.0;
			for(k = i; k < n; k++){
				sum += aa[i + k*lda] * aa[k + j*lda];
			}
		}	  
	   	  aa[i + j*lda] = sum;
	  }
  }  

}


/* =================================================
   MAIN  
   ================================================= */

int
main (void){

  pmatrix A, Ainvers, T;		/* Matrices */
  double norm;				/* Norm */
  int n = 4;				/* Problem dimension */
  /* Chose a problem */
  A = new_4x4_matrix();
  Ainvers= new_4x4_matrix();
  T = new_zero_matrix(n, n);
  print_matrix(Ainvers);
  
  /* LR - decomposition */
  printf("Start decomposition \n");
  lr_decomp(Ainvers);
  print_matrix(Ainvers);
  
  /* Invert L and R */
  printf("Start inversion of L and R \n");
  lr_invert(Ainvers);
  print_matrix(Ainvers);

  /* Multiply  */
  printf("Multiplication of R^{-1} and L^{-1} \n");
  lr_mm(Ainvers);
  print_matrix(Ainvers);
  
  /* Test invers  */
  mm(A, Ainvers, T); 
  norm = normfrob_diff_matrix(T);
  printf("At least the test\n");
  printf("|| A - A^{-1} - Id ||_{F} = %e\n", norm); 
  print_matrix(T);
  
  /* cleaning up */
  del_matrix(A);
  del_matrix(Ainvers);
  del_matrix(T);
  
  return EXIT_SUCCESS;
}
