/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 7 - Approximation von Funktionen (Interpolation) 	 */
/* ------------------------------------------------------------- */
/*	Autoren: Christina Boerst				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/

		
  #include <stdio.h>		
  #include <stdlib.h> 
  #include <math.h>
  #include <assert.h>
  #include <GL/freeglut.h>
  #include <GL/glut.h>
  
  #include "interpolation.h"
  
  #ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif	

  
/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
	
	pinterpolation interA;
	pinterpolation interB;
	pinterpolation interC;
	function f;
	void* data;
	double p;
	int width_C;

  
/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
double exampleExp(double x, void* data){
	double lambda=((double*)data)[0];
	return lambda*exp(-x*x/2);
}

double examplePoly(double x, void* data){
	return 1/(1+x*x);
}


void mix(pinterpolation inter, double p, pinterpolation A, pinterpolation B){
	
	int m = inter->m;
	double* xA = A->xi;
	double* xB = B->xi;
	double* xi = inter->xi;
	
	for(int i=0; i<m+1; i++){
		xi[i]=p*xA[i]+(1.0-p)*xB[i];
	}
	
}

/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/

static void
reshape_A(int width, int height){
	
    glViewport(0, 0, width, height); 
	
	glLoadIdentity();
	if(width > height) glScalef((double) height/width, 1.0, 1.0);
    else glScalef(1.0, (double) width/height, 1.0);
}


static void
reshape_B(int width, int height){
	
    glViewport(0, 0, width, height); 
	
	glLoadIdentity();
	if(width > height) glScalef((double) height/width, 1.0, 1.0);
    else glScalef(1.0, (double) width/height, 1.0);
}

static void
reshape_C(int width, int height){
	
    glViewport(0, 0, width, height); 
	
	width_C=width;
	
	glLoadIdentity();
	if(width > height) glScalef((double) height/width, 1.0, 1.0);
    else glScalef(1.0, (double) width/height, 1.0);
}

static void displayFunc(function f, double xmin, double xmax, double ymin, double ymax){
	
	glBegin(GL_LINE_STRIP);
	
	for(double x=-1.0; x<=1.0; x+=0.01){
		glVertex2d(x,(f((x/2+0.5)*(xmax-xmin)+xmin,data)-ymin)/(ymax-ymin)*2-1);
	}
	
	glEnd();
	
}

static void displayInterpolation(pinterpolation inter,double xmin, double xmax, double ymin, double ymax){
	
	/* Interpolation graph */
	glBegin(GL_LINE_STRIP);
	
	for(double x=-1.0; x<=1.0; x+=0.01){
		glVertex2d(x,(eval_interpolation_polynomial(inter,(x/2+0.5)*(xmax-xmin)+xmin)-ymin)/(ymax-ymin)*2-1);
	}
	
	glEnd();
	
	/* nodes */
	
	int m = inter->m;
	double* x = inter->xi;
	double* y = inter->f;
	double dx=0.01;
	
	for(int i=0; i<m+1; i++){
		glBegin(GL_LINE_LOOP);
		glVertex2d( dx+(x[i]-xmin)/(xmax-xmin)*2-1,    (y[i]-ymin)/(ymax-ymin)*2-1);
		glVertex2d(    (x[i]-xmin)/(xmax-xmin)*2-1, dx+(y[i]-ymin)/(ymax-ymin)*2-1);
		glVertex2d(-dx+(x[i]-xmin)/(xmax-xmin)*2-1,    (y[i]-ymin)/(ymax-ymin)*2-1);
		glVertex2d(    (x[i]-xmin)/(xmax-xmin)*2-1,-dx+(y[i]-ymin)/(ymax-ymin)*2-1);
		glEnd();
	}
	
}

static void
display_A(){
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	double xmin=-10;
	double xmax=10;
	double ymin=-0.2;
	double ymax=1.2;
	
	/* original function */
	glColor3f(0.2,0.2,0.2);
	displayFunc(f,xmin,xmax,ymin,ymax);
	
	/* polynomial function */
	glColor3f(0.9,0.0,0.0);
	displayInterpolation(interA,xmin,xmax,ymin,ymax);
	
	
	glFlush();
    glutSwapBuffers();
}

static void
display_B(){
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	double xmin=-10;
	double xmax=10;
	double ymin=-0.2;
	double ymax=1.2;
	
	/* original function */
	glColor3f(0.2,0.2,0.2);
	displayFunc(f,xmin,xmax,ymin,ymax);
	
	/* polynomial function */
	glColor3f(0.0,0.7,0.0);
	displayInterpolation(interB,xmin,xmax,ymin,ymax);
	
	glFlush();
    glutSwapBuffers();
}

static void
display_C(){
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	double xmin=-10;
	double xmax=10;
	double ymin=-0.2;
	double ymax=1.2;
	
	/* original function */
	glColor3f(0.2,0.2,0.2);
	displayFunc(f,xmin,xmax,ymin,ymax);
	
	/* equidistant interpolation */
	glColor3f(0.9,0.0,0.0);
	displayInterpolation(interA,xmin,xmax,ymin,ymax);
	
	/* tschernobyl interpolation*/
	glColor3f(0.0,0.7,0.0);
	displayInterpolation(interB,xmin,xmax,ymin,ymax);
	
	/* tschernobyl interpolation*/
	glColor3f(0.9,0.7,0.0);
	displayInterpolation(interC,xmin,xmax,ymin,ymax);
	
	glFlush();
    glutSwapBuffers();
}

void motion_C(int position_x, int position_y){
	p=(double)position_x/(double)width_C;
	if(p<0.0) p=0.0;
	if(p>1.0) p=1.0;
	mix(interC, p, interA, interB);
	eval_interpolated_values(interC, f, data);
	newton_divided_differences(interC);
	
	glutPostRedisplay();
}

int
main(int argc,char **argv){
	
	double lambda=1.0;
	data=(void*)&lambda;
	int m=50;
	double a=-5;
	double b=5;
	f=examplePoly;
	
	interA=new_interpolation(m);
	setup_aequidistant_interpolationpoints(interA, a, b);
	eval_interpolated_values(interA, f, data);
	newton_divided_differences(interA);
	
	interB=new_interpolation(m);
	setup_chebyshev_interpolationpoints(interB, a, b);
	eval_interpolated_values(interB, f, data);
	newton_divided_differences(interB);
	
	p=0.5;
	interC=new_interpolation(m);
	mix(interC, p, interA, interB);
	eval_interpolated_values(interC, f, data);
	newton_divided_differences(interC);
	
	
	
	
	
    /*-------------------------------------------------
	*	 Initalize GUI
	*--------------------------------------------------*/
    glutInit(&argc, argv);
	
	/* Window A */
	glutCreateWindow("Equidistant");
	
	glutPositionWindow(100, 100);
	glutReshapeWindow(450, 450);
	glClearColor(0.3f, 0.3f, 0.3f, 0.0f);
	
	glutReshapeFunc(reshape_A);
	glutDisplayFunc(display_A);
	
	
	/* Window B */
	glutCreateWindow("Chebyshev");
	
	glutPositionWindow(600, 100);
	glutReshapeWindow(450, 450);
	glClearColor(0.3f, 0.3f, 0.3f, 0.0f);
	
	glutReshapeFunc(reshape_B);
	glutDisplayFunc(display_B);
	
	/* Window C */
	glutCreateWindow("Comparison");
	
	glutPositionWindow(1100, 100);
	glutReshapeWindow(450, 450);
	glClearColor(0.3f, 0.3f, 0.3f, 0.0f);
	
	glutReshapeFunc(reshape_C);
	glutDisplayFunc(display_C);
    glutMotionFunc(motion_C);
	
	glutMainLoop();
	
  
  return EXIT_SUCCESS;
}