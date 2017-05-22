/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 5 - Simulation Loesung der Wellengleichung 	 */
/* ------------------------------------------------------------- */
/*	Autoren: Christina Boerst				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

#ifndef GRIDFUNC1D_H
#define GRIDFUNC1D_H

#include "basic.h"

/*
  Create two structs (with pointer and pointer to const struct ) :
  -  "grid1d" (the grid)
  -  "gridfunc1d" (functions on grid)
*/

typedef struct _grid1d grid1d;
typedef grid1d *pgrid1d;
typedef const grid1d *pcgrid1d;

typedef struct _gridfunc1d gridfunc1d;
typedef gridfunc1d *pgridfunc1d;
typedef const gridfunc1d *pcgridfunc1d;


struct _grid1d{
  unsigned int n;
  double h;
};

struct _gridfunc1d{
  pcgrid1d g; 
  unsigned int d;
  double *x;
};



/*-------------------------------------------------
*	 Create and destroy a grid 
*--------------------------------------------------*/

pgrid1d
new_grid1d(unsigned int n);

void
del_grid1d(pgrid1d g);


/*------------------------------------------------
*      Create and destroy a grid function
*-------------------------------------------------*/

pgridfunc1d
new_gridfunc1d(pcgrid1d g);

void
del_gridfunc1d(pgridfunc1d u_h);


/*------------------------------------------------
*	    Utility functions
*-------------------------------------------------*/

/* Set Zero */
void
zero_gridfunc1d(pgridfunc1d u_h);

/* Sinus wave for a parts of boundary */

void
left_boundary_gridfunc1d(pgridfunc1d u_h, double t);

void
right_boundary_gridfunc1d(pgridfunc1d u_h, double t);


  #endif






