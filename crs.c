/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 		Serie 10 - Schwachbesetzte Matrizen 		 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker			 		 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "crs.h"
#include "miniblas.h"


 /* ------------------------------------------------------------
 * Constructor and destructor
 * ------------------------------------------------------------ */

pcrsmatrix
new_crsmatrix( unsigned int rows, unsigned int cols){

  pcrsmatrix crs;
  
  crs = (pcrsmatrix) calloc(1, sizeof(crsmatrix));
  crs->row = rows;
  crs->col = cols;
  
  crs->Aa = NULL;
  crs->Ai = NULL;
  crs->Aj = NULL;

  return crs;
}

void
del_crsmatrix(pcrsmatrix crs){


  if(crs->Aa){
    free(crs->Aa);
  }
  if(crs->Ai){
    free(crs->Ai);
  }
  if(crs->Aj){
//    free(crs->Aj);
  }
  
  free(crs);
}

 /* ------------------------------------------------------------
 * Set up example matrix and right hand side
 * ------------------------------------------------------------ */

pcrsmatrix
setup_poisson(unsigned int m){
    // TODO
    // überarbeiten für m > 1 
    assert(m > 1);
    double hinsq = (m+1) * (m+1); // = h ^ (-2)
    pcrsmatrix A = new_crsmatrix(m, m);

    unsigned int *Ai = calloc(m+1, sizeof(unsigned int));
    Ai[0] = 0;
    Ai[1] = 2;
    for(int i = 2; i < m; i++){
        Ai[i] = Ai[i-1] + 3;
    }
    Ai[m] = Ai[m-1] + 2;
    printf("%d\n",Ai[m]);
    double *Aa = calloc(Ai[m], sizeof(double));
    Aa[0] = 2.0 * hinsq;
    Aa[1] = -1.0 * hinsq;
    for(int i = 2; i < Ai[m]-2; i += 3){   //Ai[m]
        Aa[i] = -1.0 * hinsq;
        Aa[i+1] = 2.0 * hinsq;
        Aa[i+2] = -1.0 * hinsq;
    }
    Aa[Ai[m]-2] = -1.0 * hinsq;
    Aa[Ai[m]-1] = 2.0 * hinsq;
    
    
    for(int i = 0; i < Ai[m]; i++){
        printf("Aa[%d] = %f\n",i,Aa[i]);
    }
    
    // überarbeiten für m > 1
    unsigned int *Aj = calloc(Ai[m], sizeof(unsigned int));
    Aj[0] = 0;
    Aj[1] = 1;    
    Aj[2] = 0;
//     Aj[3] = 1;
//     Aj[4] = 2;
    
    for(int i = 3; i < Ai[m]; i++){
        Aj[i] = Aj[i-3] + 1;
    }
    
    //check
    for(int i = 0; i < Ai[m]; i++){
        printf("Aj[%d] = %d\n",i,Aj[i]);
    }
    
    // regular:
    A->Aa = Aa;
    A->Ai = Ai;
    A->Aj = Aj;
    
    return A;
}


void 
set_righthandside(int m, function f, void *data, pvector b){

  double h = 1.0/(m+1);
  unsigned int i;

  for(i = 0; i < b->rows; i++){
    b->x[i] = f((i+1)*h, data);

  }

}


 /* ------------------------------------------------------------
 * Auxiliary functions
 * ------------------------------------------------------------ */

void
print_crs(pcrsmatrix crs){

  unsigned int nnze = crs->Ai[crs->row];
  unsigned int i;


  if(nnze >1){
    printf("Non-zero entries\n");
    printf("(%f,  ", crs->Aa[0]);
    for(i = 1; i < nnze-1; i++){
      printf("%f,  ", crs->Aa[i]);
    }
    printf("%f)\n", crs->Aa[nnze-1]);
    printf("\n");
    printf("Entries per row\n(%d", crs->Ai[0]);
    for(i = 1; i <= crs->row; i++){
      printf(", %d", crs->Ai[i]);
    }
    printf(")\n");
    printf("Column for entries\n");
    printf("(%d, ", crs->Aj[0]);
    for(i = 1; i < nnze-1; i++){
      printf("%d, ", crs->Aj[i]);
    }
    printf("%d)\n", crs->Aj[nnze-1]);

  }
  else{
    printf("Matrix\n");
    printf("(%f ) \n", crs->Aa[0]);
  }


}


void
mvm_crs(pcrsmatrix crs, pvector x, double alpha, pvector b){
    unsigned int *ai = crs->Ai;
    unsigned int *aj = crs->Aj;
    double *aa = crs->Aa;
    double *xx = x->x;
    double *bx = b->x;
    unsigned int row = crs->row;
    
    assert(crs->col == x->rows);
    assert(crs->row == b->rows);
    
    for(int i = 0; i < row; i++){
        bx[i] *= alpha;
        for(int j = ai[i]; j < ai[i+1]; j++){
            bx[i] += xx[aj[j]] * aa[j]; 
        }
    }
}

void
richardson_iteration(pcrsmatrix crs, pvector x, double theta, pvector b, double eps){
    unsigned int dim = crs->row;
    pvector h = new_vector(dim);
    do{
        for(int i = 0; i < dim; i++){
            h->x[i] = b->x[i];
        }
          
        mvm_crs(crs, x, -1.0, h);
        axpy(dim, -theta, h->x, 1.0, x->x, 1.0);
    }while(nrm2(dim,h->x,1.0) > eps);

}

