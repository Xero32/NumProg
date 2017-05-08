
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "matrix.h"
#include "miniblas.h"

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

pmatrix
new_matrix(int rows, int cols){

  pmatrix a;

  a = (pmatrix) malloc(sizeof(matrix));
  a->a = (double *) malloc(sizeof(double) * rows * cols);
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;
  
  return a;
}

pmatrix
new_zero_matrix(int rows, int cols){

  pmatrix a;

  a = calloc(1, sizeof(matrix));
  a->a = calloc(rows * cols, sizeof(double));
  a->rows = rows;
  a->cols = cols;
  a->ld = rows;

  return a;
}

pmatrix
new_identity_matrix(int rows){

  pmatrix a;
  int i;

  a = new_zero_matrix(rows, rows);

  for(i = 0; i < rows; i++){
    a->a[i + i *rows] = 1.0;
  }

  return a;
}

void
del_matrix(pmatrix a){

  free(a->a);
  free(a);
}


pvector
new_vector(int rows){

  pvector x;

  x = (pvector) malloc(sizeof(vector));
  x->x = (double *) malloc(sizeof(double) * rows);
  x->rows = rows;

  return x;
}

pvector
new_zero_vector(int rows){
	
  pvector x;
  int i;
 
  x = new_vector(rows);
  
  for(i = 0; i < rows; i++){
	x->x[i] = 0.0;  
  }
	
  return x;	
}

void
del_vector(pvector x)
{
  free(x->x);
  free(x);
}


pvector
matrix_col(pmatrix a, int i){
  
  pvector x;
  
  assert(i < a->rows);
  x =  new_vector(a->rows);
  x->x =a->a +(i*(a->rows)); 

  return x;
}


/* Submatrices 				*/

pmatrix //TODO
init_sub_matrix(pmatrix asub, pmatrix a, int rows, int roff, int cols, int coff){
    
//     asub->a = a->a+roff+coff*a->ld;
//     asub->ld=a->ld;
//     asub->rows = rows;
//     asub->cols = cols;
//     printf("matrix asub:\n");
//     print_matrix(asub);
//     
//     printf("matrix a:\n");
//     print_matrix(a);
    double *A = a->a;
    int lda = a->ld;
    double *AS = asub->a;
    for(int j = coff; j < cols+coff; j++){
        for(int i = roff; i < rows+roff; i++){
            int k = i - roff;
            int m = j - coff;
            //asub->a+k+m*rows = a->a+i+j*lda;
            AS[k+m*rows] = A[i+j*lda];
        }
    }
//     printf("matrix asub:\n");
//     print_matrix(asub);
    return asub;
}

/* ------------------------------------------------------------
 * Example matrix
 * ------------------------------------------------------------ */

pmatrix 
new_diaghilbert_matrix(int rows){

  pmatrix a;
  double *aa;
  double sum;
  int lda;
  int i, j;

  a = new_matrix(rows, rows);
  aa = a->a;
  lda = a->ld;

  for(j=0; j<rows; j++) {
    sum = 1.0;
    for(i=0; i<j; i++) {
      aa[i+j*lda] = 1.0 / (1.0 + i + j);
      sum += fabs(aa[i+j*lda]);
    }
    for(i=j+1; i<rows; i++) {
      aa[i+j*lda] = 1.0 / (1.0 + i + j);
      sum += fabs(aa[i+j*lda]);
    }
    aa[j+j*lda] = sum;
  }

  return a;
}

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

/* ------------------------------------------------------------
 * Simple utility functions
 * ------------------------------------------------------------ */

void
clear_matrix(pmatrix a){

  double *aa = a->a;
  int rows = a->rows;
  int cols = a->cols;
  int lda = a->ld;
  int i, j;

  for(j=0; j<cols; j++)
    for(i=0; i<rows; i++)
      aa[i+j*lda] = 0.0;
}

void
clear_vector(pvector x){

  double *xx = x->x;
  int rows = x->rows;
  int i;
  
  for(i=0; i<rows; i++)
    xx[i] = 0.0;
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

