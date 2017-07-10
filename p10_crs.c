/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/*   	Serie 10 - Schwachbesetzte Matrizen (Richardson)	 */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker	 				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

  #include <stdio.h>		
  #include <stdlib.h> 
  #include <math.h>
  #include "crs.h"

  #ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif

/* right hand side */
static double 
sinus(double x, void *data){
	
  double *dat = (double *) data;
  double fx;	
  double mu = dat[0];
	  
  fx = sin(M_PI * mu * x);	
  
  return fx;
}

void
WriteToFile(FILE *F, int M, double Theta, double Diff){
    if(F) fprintf(F, "%d\t%f\t%e\n", M, Theta, Diff);    
}

int 
main(void){
  
  pcrsmatrix crs;
  unsigned int m;
  pvector x, b;
  unsigned int i;
  double data[1];
  double theta, h, lambda;
  
//   int ctr = 0;
  double diff;
  m = 2;					/* problem dimension */
  data[0] = 1.0;
  theta = 0.01;
  
    
//   double theta0 = 0.01;
//   FILE *f;
//   f= fopen("richardson7.txt", "a");
//   do{
//     switch(ctr){
//       default: break;
//       case 0: ctr++; theta = theta0; break;
//       case 1: ctr++; theta = theta0 * 2; break;
//       case 2: ctr++; theta = theta0 * 4; break;
//       case 3: ctr++; theta = theta0 * 6; break;
//       case 4: ctr++; theta = theta0 * 8; break;
//       case 5: ctr++; theta = theta0 *= 10; break;
//       case 6: ctr = 1; 
//             if(theta0 >= 0.02){
//                 theta0 = 0.0000001;
//                 m++;
//             }
//             continue;
//     }
  
  h = 1.0/(m+1);
  printf("Solve poisson equation with Richardson iteration\n");
  printf("-----------------------------------------\n");
 
  /* set up problem */  
  crs = setup_poisson(m);  	/* matrix */
  b = new_zero_vector(m);	/* right hand side */
  x = new_zero_vector(m);	/* solution */
  set_righthandside(m, sinus, data, b);
  
  printf("Right hand side \n");
    for(i = 0; i < m; i++){
    printf("%f\n", b->x[i]);
  }
  /* solve with Richardson iteration */
  richardson_iteration(crs, x, theta, b, 1.e-7);
  printf("-----------------------------------------\n");
  printf("Solution\n");
  	lambda = sinus(h *0.5, data);		/* sin(k*h/2)*/
    lambda *= lambda;					/* sin^{2} */
 	lambda *= 4.0/h/h;					/* 4.0/(h*h) * sin */
  for(i = 0; i < m ; i++){
    printf("Richardson %f \t exakt %f \t difference %e \n", x->x[i],  b->x[i]/ lambda, diff = fabs(x->x[i] - b->x[i]/ lambda));
  }
 
//  WriteToFile(f, m, theta, diff);
 

 
//  printf("print poisson matrix:\n");
//  print_crs(crs);
  /* Cleaning up */
  del_crsmatrix(crs);
  del_vector(x);
  del_vector(b);
  
//   }while(m < 100);
//   if(f) fclose(f);
  return EXIT_SUCCESS;
}
