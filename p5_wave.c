/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 5 - Simulation Loesung der Wellengleichung 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 						 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

		
#include <stdio.h>		
#include <stdlib.h> 
#include <GL/freeglut.h>
#include <GL/glut.h>
  
#include "gridfunc1d.h"
#include "leapfrog1d.h"	
  
/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/

/* Some nice ideas for global values, but feel free to change them */

/*	Wave equation and grid		*/							
pgridfunc1d u[2], v[2];  
unsigned int current = 0;		/* switch between grid functions */
double data[2];				/* data 'c' and left or right wave */
double t = 1.25;			/* time used to create a start wave */
double delta;				/* incremenet */
unsigned int step;			/* to find a good relation between increment size (therefore accuracy) 
 					   update rate for glut. */
  
 
/* reshape function (simple 2D without frills!) */ 
static void
reshape_wave(int width, int height){
  
  /* ---------------------------------------------- */ 
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */ 
  /*     T       O     O     D     D     O     O    */ 
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */ 
  /* ---------------------------------------------- */ 	

}
 
 
/* maybe you want to use different colors ? ;-) */
static void
color(double coefficient){
	  

	  
} 
  
  
/* and now... the content =)
   draw the string and the deflections, maybe with some nice colors
   and think about scaling !*/

static void
display_wave(){
	  
  /* ---------------------------------------------- */ 
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */ 
  /*     T       O     O     D     D     O     O    */ 
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */ 
  /* ---------------------------------------------- */   
  
}
  
/* start new waves, leave and so on.... there are a lot of possibilities 
   Think about a way to change the start point for a new wave (left or right) 
   and how to reset the simulation */   
static void
key_wave(unsigned char key, int x, int y){

  /* ---------------------------------------------- */ 
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */ 
  /*     T       O     O     D     D     O     O    */ 
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */ 
  /* ---------------------------------------------- */  
	  
}

  
/* it's necessary to compute new values and redisplay them after a while....*/
static void
timer_wave(int val){
	  
  /* ---------------------------------------------- */ 
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */ 
  /*     T       O     O     D     D     O     O    */ 
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */ 
  /* ---------------------------------------------- */ 	  
	  
}
  
  
/* last but not least the main function 
   don't forget to set 'm' 'c' and the number
   of points 'n' */   
int
main(int argc,char **argv){
  
   /* ---------------------------------------------- */ 
  /*                                                */
  /* T T T T T     O O       D D           O O      */
  /*     T        O   O      D   D        O   O     */
  /*     T       O     O     D     D     O     O    */ 
  /*     T       O     O     D     D     O     O    */ 
  /*     T        O   O      D   D        O   O     */
  /*     T         O O       D D           O O      */
  /*                                                */ 
  /* ---------------------------------------------- */ 
  
  return EXIT_SUCCESS;
}
