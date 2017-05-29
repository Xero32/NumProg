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
double t = 2.0;			/* time used to create a start wave */
double delta = 0.01;				/* incremenet */
// smaller?
unsigned int step;			/* to find a good relation between increment size (therefore accuracy) 
 					   update rate for glut. */
unsigned int z = 1;

// delta = 0.05;
/* reshape function (simple 2D without frills!) */ 
static void
reshape_wave(int width, int height){
  
    glViewport(0, 0, width, height);
    windowWidth=width;
    windowHeight=height;
//     glMatrixMode(GL_PROJECTION);  // uncertain about use
    glLoadIdentity();
//     if(width > height){								
//        glScalef((double) height/width, 1.0, 1.0);
//     }
//     else{												
//        glScalef(1.0, (double) width/height, 1.0);
//     }	
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
    glClearColor(0.3f,0.3f,0.3f,0.0f);
    if(z){
        glPushMatrix();
        glLoadIdentity();
        glBegin(GL_LINES);
        color(1.0);
            glVertex2f(-1.0,0.0);
            color(1.0);
            glVertex2f(1.0,0.0);
            color(1.0);
        glEnd();
        glFlush();
        glutSwapBuffers();
    }
    
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
            s = (fabs(s) < fabs(y[i])) ? fabs(y[i]) : s;
        }
        
	glBegin(GL_LINE_STRIP);
	
	  
	for(int i=0; i<n; i++){
            color(0.0);
            glVertex2f((double)i/(double)n * 2.0 - 1, 0.6*y[i]/s);
	}
	
	glEnd();
// 	*/
	
	pgridfunc1d pos=u[current];
	y=pos->x;
	n=pos->d;
	
	glBegin(GL_LINE_STRIP);
	
        s = 0.0;
        for(int i = 0; i < n; i++){
            s = (fabs(s) < fabs(y[i])) ? fabs(y[i]) : s;
        }
	 
	for(int i=0; i<n; i++){
            color(1.0); 
            glVertex2f((double)i/(double)n * 2.0 - 1, 0.6*y[i]/s);
	}
    
    glEnd();
    glFlush();
    glutSwapBuffers();
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
    if(step == 3){
        glutPostRedisplay();	
        step = 0;
    }
    glutTimerFunc(15, timer_wave, 0);  

    t += delta;
    step++;
//     printf("t: %f,\tz: %f\n",t,z);
    
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
          case 'l':             //left_boundary   
                                data[0] = 'l';
                                t = 0.0;
                                z = 0;
                                printf("LEFT\n");
                                break;
          case 'r':             //right_boundary
                                data[0] = 'r';
                                t = 0.0;
                                z = 0;
                                printf("RIGHT\n");
                                break;
          case 'z':             zero_gridfunc1d(u[current]);
                                zero_gridfunc1d(v[current]);
                                t = 2.0;
                                z = 1;
                                
                                
                                glPushMatrix();
                                glLoadIdentity();
                                glBegin(GL_LINES);
                                    color(1.0);
                                    glVertex2f(-1.0,0.0);
                                    color(1.0);
                                    glVertex2f(1.0,0.0);
                                color(1.0);
                                glEnd();
                                printf("ZERO\n");
                                glFlush();
                                glutSwapBuffers();

                                break;
	  default :
				break;	 		
	  }
	  
}

  

  
/* last but not least the main function 
   don't forget to set 'm' 'c' and the number
   of points 'n' */   

double
f(double k, double delta){
    return 0.065 * exp(-0.055 * k) + 0.007 - delta;
}



int
main(int argc,char **argv){
    data[1] = 0.15;
    
  if(argc == 3){
    data[1] = atof(argv[1]);//c
    printf("c: %f\n",data[1]);
    delta = atof(argv[2]);
    printf("delta: %f\n",delta);
  }else if(argc == 2){
//       printf("CHECK\n");
      data[0] = argv[1][0];
//       printf("data: %f\n",data[0]);
    if(data[0] == 0x68){
        printf("Usage info:\nSet perturbation on either left or right end by pressing:\n\t 'l': left\n\t 'r': right\n\n");
        printf("You can reset the wave by pressing 'z'\nResetting may not work properly at first try.\nIn that case press 'z' repeatedly.\n\n");
        printf("Change Hooke's constant 'c' with input value.\n");
        printf("Change time increment 'delta' with second input value.\n\n");
        printf("For c and delta there is an exponential relation, which cannot be exceeded for the program to work properly.\n\n");
        
        exit(EXIT_SUCCESS);
    }else{
        data[1] = atof(argv[1]);
        printf("c: %f\n",data[1]);
    }
  }else{
        data[0] = 'l';
        data[1] = 0.1;
        delta = 0.007;
  }
    
    int n = 300;
    double k = 2.0 * data[1] * data[1] * (n+1.0) * (n+1.0) * delta;
    printf("f: %f\n",f(k,delta));
    if( f(k,delta) <  0 ){
            printf("Wrong ratio of delta and c chosen\nChoose lower values\n");
            exit(EXIT_SUCCESS);
    }
    printf("Actual values:\n\tc: %f\tlr: %f\tdelta: %f\n",data[1],data[0],delta);
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
  zero_gridfunc1d(u[current]);
  glutKeyboardFunc(key_wave);
  glutTimerFunc(50, timer_wave, 0);
  glutMainLoop();
  return EXIT_SUCCESS;
}
