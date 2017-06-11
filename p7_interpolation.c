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
int m = 10;
// int n = 50;
int n;
// const double lambda = 1.0;
double data[2] = {1.0,1.0};

pinterpolation inter0, inter1;
/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
double
f(double x, double lambda){
    assert(lambda);
    return lambda * exp(-0.5 * pow(x,2));
}

double
g(double x){
    assert(1 + pow(x,2)); 
    return 1 / (1 + pow(x,2));
}
/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/
void
display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glPushMatrix();
    
    
    // draw function g
    inter0 = new_interpolation(m);
    double *y = inter0->f;
    double *xi = inter0->xi;
    setup_aequidistant_interpolationpoints(inter0, -1.0, 1.0);
    eval_interpolated_values(inter0, (function)g, data); // ho to use data?
// draw interpolation
    glBegin(GL_LINE_STRIP);
    glColor3f(0.9,0.0,0.0);
    for(int i = 0; i < m+1; i++){
        glVertex2f(xi[i], y[i]-0.25);
    }
    glEnd();
    
    glBegin(GL_LINE_STRIP);
    glColor3f(0.5,0.5,0.5);
//  draw proper function
    for(double i = -1.0; i < 1.0; i += 0.01){
        glVertex2f(i, 1 / (1 + pow(i,2)) - 0.25 );
    }
    glEnd();    
    
    // draw function f
    inter1 = new_interpolation(m);
    double *y1 = inter1->f;
    double *x1 = inter1->xi;
    setup_aequidistant_interpolationpoints(inter1, -1.0, 1.0);
    eval_interpolated_values(inter1, (function)f, data);
    //draw interpolation
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.8,0.0);
    for(int i = 0; i < m+1; i++){
        glVertex2f(x1[i], y1[i] - 0.25);
        printf("i: %d, xi[i]: %f, y[i]: %f\n",i,x1[i],y1[i]);
    }
    glEnd();
    
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.0);
    //draw proper function
    for(double i = -1.0; i < 1.0; i += 0.01){
        glVertex2f(i, 1.0 * exp(-0.5 * pow(i,2)) - 0.25);
    }
    glEnd();
    
    
    glFlush();
    glutSwapBuffers();
}

void
reshape(int width, int height){
    glViewport(0, 0, width, height); 
    glLoadIdentity();
    if(width > height){								
        glScalef((double) height/width, 1.0, 1.0);
    }
    else{												
        glScalef(1.0, (double) width/height, 1.0);
    }
}

/* If necessary use the pictures as orientation,
   but feel free to create something different =)*/

void
key(unsigned char key, int x, int y){
    switch(key){
        case 27: exit(EXIT_SUCCESS);
    }
}

int
main(int argc, char **argv){
    assert(n && m);
    glutInit(&argc, argv);
    glutCreateWindow("Interpolation");
    glutPositionWindow(150, 100);
    glutReshapeWindow(800, 800);  
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutMainLoop();
    return EXIT_SUCCESS;
}