/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 		Serie 10 - Schwachbesetzte Matrizen 		 */
/* ------------------------------------------------------------- */
/*	Autoren: 				 		 */
/*	Versionsnummer:						 */
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
 //   free(crs->Aj);
  }
  
  free(crs);
}

 /* ------------------------------------------------------------
 * Set up example matrix and right hand side
 * ------------------------------------------------------------ */

pcrsmatrix
setup_poisson(unsigned int m){
 
	pcrsmatrix crs = new_crsmatrix(m, m);  
	double h = (m+1)*(m+1);
	unsigned int* Ai = calloc(m+1, sizeof(unsigned int));
	Ai[0]=0;
	Ai[1]=2;
	for(int i=2;i<m;i++){
		Ai[i] = Ai[i-1]+3;
	}
	Ai[m] = Ai[m-1] +2;
	
	double* Aa = calloc(m*3-2, sizeof(double));
	unsigned int* Aj = calloc(m*3-2, sizeof(unsigned int));
	
	Aa[0]=2*h;
	Aa[1]=-1*h;
	for(int i=2;i<m*3-4;i+=3){
		Aa[i]=-1*h;
		Aa[i+1]=2*h;
		Aa[i+2]=-1*h;
	}
	Aa[m*3-4]=-1*h;
	Aa[m*3-3]=2*h;
	
	Aj[0]=0;
	Aj[1]=1;
	
	for(int i=1;i<m-1;i++){
		Aj[3*i-1]=i-1;
		Aj[3*i]=i;
		Aj[3*i+1]=i+1;
	}
	
	Aj[3*m-4]=m-2;
	Aj[3*m-3]=m-1;
	
	crs->Ai =Ai;
	crs->Aj =Aj;
	crs->Aa =Aa;
	
	return crs;
	
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
	
	unsigned int cr = crs->row;
	unsigned int cc = crs->col;
	unsigned int* ci = crs->Ai;
	unsigned int* cj = crs->Aj;
	double* ca = crs->Aa;
	double* bx = b->x;
	double* xx = x->x;

	int xr = x->rows;
	int br = b->rows;
	
	unsigned int i;
	unsigned int j;
	
	assert(cc==xr);
	assert(cr==br);
	
	for(i=0;i<cr;i++){
		bx[i] *= alpha;
		for(j=ci[i];j<ci[i+1];j++){
			bx[i] += ca[j]*xx[cj[j]];
		}
	}
}

void
richardson_iteration(pcrsmatrix crs, pvector x, double theta, pvector b, double eps){
	
	unsigned int rows = crs->row;
	pvector h = new_vector(rows);
	int i;
	
	
	for(i=0;i<rows;i++){
		h->x[i] = b->x[i];
	}
	mvm_crs(crs, x, -1.0, h);
	axpy(rows, -theta, h->x, 1, x->x, 1);
	
	while(nrm2(rows,h->x,1)>eps){
	for(i=0;i<rows;i++){
		h->x[i] = b->x[i];
	}
	mvm_crs(crs, x, -1.0, h);
	axpy(rows, -theta, h->x, 1, x->x, 1);
	}
	//for(i=0;i<rows;i++){
	//	printf("%f\n",h->x[i]);
	//}
  
}

