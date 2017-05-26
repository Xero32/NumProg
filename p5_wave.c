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
int n;
pgridfunc1d u[2], v[2]; 
int windowWidth, windowHeight;
unsigned int current = 0;		/* switch between grid functions */
double data[2];				/* data 'c' and left or right wave */
//data[1] = c;
double t = 0.0;			/* time used to create a start wave */
double delta = 0.05;				/* incremenet */
unsigned int step;			/* to find a good relation between increment size (therefore accuracy) 
 					   update rate for glut. */
// delta = 0.05;
/* reshape function (simple 2D without frills!) */ 
static void
reshape_wave(int width, int height){
  
    glViewport(0, 0, width, height);
    windowWidth=width;
    windowHeight=height;
//     glMatrixMode(GL_PROJECTION);  // uncertain about use
    glLoadIdentity();
    if(width > height){								
       glScalef((double) height/width, 1.0, 1.0);
    }
    else{												
       glScalef(1.0, (double) width/height, 1.0);
    }	
}
 
 
/* maybe you want to use different colors ? ;-) */
static void //possibly ignore
color(double coefficient){
    
	glColor3f(coefficient,coefficient,coefficient);
	  
} 
  
  
/* and now... the content =)
   draw the string and the deflections, maybe with some nice colors
   and think about scaling !*/

static void
display_wave(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
//     glClearColor(0.0,0.0,0.0,1.0);
//     glPushMatrix();
//     glBegin(GL_LINE_STRIP); 
//     color(0.2);
//     glColor3f(0.0,0.1,0.2);
    
//     double *y = u[current]->x;
//     glutSolidTeapot(0.5);
//     for(double j = 0; j <= 1.0; j+=0.07){
//         glVertex2f(-0.5,j);
//         glVertex2f(0.5,j);
//         
//         glVertex2f(-0.5,j);
//         glVertex2f(0.5,-j);
// }
//        glVertex2f((double)2*(i/n) - 1, y[i]);     //scale!!!
//        glColor3f(0.8,0.1,0.2);
//     glEnd();
    
    //alternatively, Tim's for-loop:
//     for(int i=0; i<=1.0; i+=0.07){
// 	  glVertex2f(-0.9,i);
// 	  glVertex2f(0.9,i);
// 	  
// 	  glVertex2f(-0.9,-i);
// 	  glVertex2f(0.9,-i);
// 	}
//     
    
    
    /*-------------------------------------------------
	*	 Draw current v
	*--------------------------------------------------*/
//         int n = u[0]->d;
	
	
// 	/*
	pgridfunc1d vel=v[current];
	double* y=vel->x;
        int n=vel->d;
	
        double s = 0.0;
        for(int i = 0; i < n; i++){
            s = (s < y[i]) ? y[i] : s;
        }
        
	glBegin(GL_LINE_STRIP);
	
	  
	for(int i=0; i<n; i++){
            color(0.0);
            glVertex2f((double)i/(double)n * 2.0 - 1, 0.1*y[i]/s);
	}
	
	glEnd();
// 	*/
	
	pgridfunc1d pos=u[current];
	y=pos->x;
	n=pos->d;
	
	glBegin(GL_LINE_STRIP);
	
        s = 0.0;
        for(int i = 0; i < n; i++){
            s = (s < y[i]) ? y[i] : s;
        }
	 
	for(int i=0; i<n; i++){
            color(1.0); 
            glVertex2f((double)i/(double)n * 2.0 - 1, 0.3*y[i]/s);
	}
    
    glEnd();
    glFlush();
    glutSwapBuffers();
}
  
/* start new waves, leave and so on.... there are a lot of possibilities 
   Think about a way to change the start point for a new wave (left or right) 
   and how to reset the simulation */   
static void
key_wave(unsigned char key, int x, int y){	
	(void) x;
	(void) y;
	
	/*Beschreibung, welche Folgen das 
	Druecken einer Taste hat.
	Hier: 'esc' fuehrt zum Beenden (Fenster wird
	geschlossen)
	Alles andere hat keine Wirkung!*/	
	  switch (key){
	  case 27 :
				exit (EXIT_SUCCESS);
				break;
	  default :
				break;	 		
	  }
	  
}

  
/* it's necessary to compute new values and redisplay them after a while....*/
static void
timer_wave(int val){
    if(!current){
        step_leapfrog1d_wave(u[0], v[0], u[1], v[1], t, delta, data);
        current = 1;
    }else{
        step_leapfrog1d_wave(u[1], v[1], u[0], v[0], t, delta, data);
        current = 0;
    }
//     printf
    if(step == 5){
        glutPostRedisplay();	
        step = 0;
    }
    glutTimerFunc(15, timer_wave, 0);  

    t += delta / 10.0;
    step++;
//     printf("t: %f,\tz: %f\n",t,z);
    
}
  
/* last but not least the main function 
   don't forget to set 'm' 'c' and the number
   of points 'n' */   


int
main(int argc,char **argv){
//   if(argc == 3){
//     data[1] = *argv[1]; //c
//     data[0] = *argv[0]; //lr
//   }
    data[1] = 0.01;
    data[0] = 'l';
    int n = 200;
    
    pgrid1d grid=new_grid1d(n);
	
	u[0]=new_gridfunc1d(grid);
	zero_gridfunc1d(u[0]);
	u[1]=new_gridfunc1d(grid);
	zero_gridfunc1d(u[1]);
	
	v[0]=new_gridfunc1d(grid);
	zero_gridfunc1d(v[0]);
	v[1]=new_gridfunc1d(grid);
	zero_gridfunc1d(v[1]);
	
    
// color(0.0);
// key_wave('s',1,2);
  glutInit(&argc, argv);
  
  glutCreateWindow("Wave");
  
  glutPositionWindow(150, 100);
  
  glutReshapeWindow(600, 600);
  glClearColor(0.5f, 0.5f, 0.9f, 0.0f);	
  glutReshapeFunc(reshape_wave);
  
  glutDisplayFunc(display_wave);
//   glClearColor(0.3f,0.3f,0.3f,0.0f);
  
  glutKeyboardFunc(key_wave);
  glutTimerFunc(50, timer_wave, 0);
  glutMainLoop();
  return EXIT_SUCCESS;
}
