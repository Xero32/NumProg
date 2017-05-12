
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


/* Submatrices 				*/

pmatrix
init_sub_matrix(pmatrix asub, pmatrix a, int rows, int roff, int cols, int coff){

  assert(a->rows >= rows+roff && a->cols >= cols+coff);
	
  asub->a=a->a+roff+coff*a->ld;
  asub->ld=a->ld;
  asub->rows=rows;
  asub->cols=cols;
  return asub;
}

ptridiag 
new_tridiag(int rows){
	
  ptridiag a;

  assert(rows >= 1);
  a = calloc(1, sizeof(tridiag));
  a->rows = rows;
  a->d = calloc(rows, sizeof(double));
  a->l = NULL;
  a->u = NULL;
  if(rows > 1) {
    a->l = calloc((rows-1), sizeof(double));
    a->u = calloc((rows-1), sizeof(double));
  }

  return a;
}

void
del_tridiag(ptridiag a){
	
  free(a->u);
  free(a->l);
  free(a->d);
  free(a);
  
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

ptridiag
new_threepointstencil(int rows){
	
  ptridiag a;
  double *ad, *al, *au;
  double h = 1.0 / (rows + 1);
  double l = -1.0 / h / h;
  double d = 2.0 / h / h;
  int i;

  a = new_tridiag(rows);
  ad = a->d;
  al = a->l;
  au = a->u;

  for(i = 0; i < rows-1 ; i++){
    al[i] = l;
    au[i] = l;
    ad[i] = d;
  }
  ad[rows-1] = d;


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

void
random_vector(pvector x)
{
  double *xx = x->x;
  int rows = x->rows;
  int i;

  for(i = 0; i < rows; i++){
    xx[i] = 2.0 * rand() / RAND_MAX - 1.0;
  }
}

double
normmax_vector(const pvector x){
	
  double norm = 0.0;
  double * xx= x->x;
  int i;

  for(i = 0; i < x->rows; i++){
    norm = (norm < fabs(xx[i]) ? fabs(xx[i]) : norm);   
  }  
	
  return norm;	
}

double
normmax_diff_vector(const pvector x, const pvector y){
	
  double error = 0.0;
  int i;
  
  assert(x->rows == y->rows);
  
  for(i = 0; i < x->rows; i++){
    error = ((fabs(x->x[i] - y->x[i])) > error ? (fabs(x->x[i] - y->x[i])) : error);
  }
   
  return error;	
}

double
norm2_diff_vector(const pvector x, const pvector y){
	
  double *diff;
  double norm;
  int i;
  
  assert(x->rows == y->rows);
  diff = (double *) malloc(sizeof(double) * x->rows);
  
  for(i = 0; i < x->rows; i++){
	diff[i] = x->x[i] - y->x[i];
  }

  norm = nrm2(x->rows, diff, 1);

  return norm;
}



void
copy_tridiag(const ptridiag a, ptridiag b){
	
  const double *ad = a->d;
  const double *al = a->l;
  const double *au = a->u;
  double *bd = b->d;
  double *bl = b->l;
  double *bu = b->u;
  int rows = a->rows;
  int i;

  assert(b->rows == rows);

  for(i = 0; i < rows- 1; i++){
	bd[i] = ad[i]; 
	bl[i] = al[i];
        bu[i] = au[i];
  }
  	bd[rows-1] = ad[rows-1]; 

}

ptridiag
clone_tridiag(const ptridiag a){
	
  ptridiag b;

  b = new_tridiag(a->rows);

  copy_tridiag(a, b);

  return b;
}

/* Matrix multiplication for tridiagonal matrices */

void
mvm_tridiag(int trans, double alpha, const ptridiag a,
	    const pvector x, pvector y){

  const double *ad = a->d;
  const double *al = a->l;
  const double *au = a->u;
  const double *xx = x->x;
  double *yx = y->x;
  int rows = a->rows;
  int i;

  assert(x->rows == rows);
  assert(y->rows == rows);

  if(rows == 0)
    return;

  if(rows == 1) {
    yx[0] += alpha * ad[0] * xx[0];
    return;
  }

  if(trans) {
    yx[0] += alpha * (ad[0] * xx[0] + al[0] * xx[1]);
    for(i = 1; i < rows-1; i++){
      yx[i] += alpha * (au[i-1] * xx[i-1] + ad[i] * xx[i] + al[i] * xx[i+1]);
    }
    yx[i] += alpha * (au[i-1] * xx[i-1] + ad[i] * xx[i]);
  }
  else {
    yx[0] += alpha * (ad[0] * xx[0] + au[0] * xx[1]);
    for(i = 1; i < rows-1; i++){
      yx[i] += alpha * (al[i-1] * xx[i-1] + ad[i] * xx[i] + au[i] * xx[i+1]);
    } 
    yx[i] += alpha * (al[i-1] * xx[i-1] + ad[i] * xx[i]);
  }
}

