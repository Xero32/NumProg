/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 8 - Approximation von Integralen (Quadratur) 	 */
/* ------------------------------------------------------------- */
/*	Autoren: 						 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

		
  #include <stdio.h>		
  #include <stdlib.h> 
  #include <math.h>
  #include <assert.h>
  #include <GL/freeglut.h>
  #include <GL/glut.h>
  
  #include "quadrature.h"
  
  #ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif	

  
/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
double a, b, areamid, areatr;
double areamid_comp, areatr_comp;
int m, n;
int quadflag;
pquadrature midquad, trquad;
double data[1] = {0.0};
/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
static double
f(double x){
    return 5.0 * exp(2.0 * x) * cos(x) / (exp(M_PI)-2.0);
}

/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/
void
draw_coord(){
    glColor3f(0.0,0.0,0.0);
    glVertex2f(-1.0,-100.0);    //y-axis
    glVertex2f(-1.0,100.0);     //y-axis
    glVertex2f(100.0,0.0);      //x-axis
    glVertex2f(-100.0,0.0);     //x-axis
    glVertex2f(0.0,-.07);        //x-mark
    glVertex2f(0.0,.07);
    glVertex2f(-2.0,-.07);       //x-mark
    glVertex2f(-2.0,.07);
    glVertex2f(-0.8,1.0);       //y-mark
    glVertex2f(-1.2,1.0);  
    glVertex2f(-0.8,-1.0);      //y-mark
    glVertex2f(-1.2,-1.0);
}

void
draw_marks(double loc, double wid){
    glVertex2f(loc, wid*0.5);
    glVertex2f(loc,-wid*0.5);  
}

void
display(){
    glutSetWindow(1);
    double aa = a-1.0;
    double bb = b-1.0;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glEnable (GL_BLEND); 
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    
    // draw function
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    for(int j = 0; j < 500; j++){
        glVertex2f((double)(j/50.0)-5.0,f((double)(j/50.0)-4.0));
    }
    glEnd();
    
    // draw single midpoint quadrature area
    glColor4f(0.0,1.0,0.0,0.4);
    glBegin(GL_POLYGON);
    glVertex2f(aa,0.0);
    glVertex2f(aa,f(M_PI * 0.25));
    glVertex2f(bb,f(M_PI * 0.25));
    glVertex2f(bb,0.0);
    glEnd();
    
    //draw single trapezoid quadrature area
    glColor4f(0.0,0.0,0.0,0.4);
    glBegin(GL_POLYGON);
    glVertex2f(aa,0.0);
    glVertex2f(aa,f(a));
    glVertex2f(bb,f(b));
    glVertex2f(bb,0.0);
    glEnd();
    
    //draw coordinate system
    //shift coordinate system, because function is boring for x < 0
    glBegin(GL_LINES);
    draw_coord();
    glEnd();
    
    glLineWidth(1.5);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0);
    draw_marks(aa,.28);
    draw_marks(bb,.28);
    glEnd();
    glLineWidth(1.0);

    glFlush();
    glutSwapBuffers();
}

void
display2(){
    glutSetWindow(2);
    double aa = a-1.0;
    double bb = b-1.0;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glEnable (GL_BLEND); 
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    
     // draw function
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    for(int j = 0; j < 500; j++){
        glVertex2f((double)(j/50.0)-5.0,f((double)(j/50.0)-4.0));
    }
    glEnd();
    
    //draw coordinate system
    //shift coordinate system, because function is boring for x < 0
    glBegin(GL_LINES);
    draw_coord();
    glEnd();
    glLineWidth(1.5);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0);
    draw_marks(aa,.28);
    draw_marks(bb,.28);
    glEnd();
    glLineWidth(1.0);
    
    
//     double *xt = trquad->xq;
    double u = 1.0 / (double)n;
    double l = b-a;
    
    //draw composite midpoint quadrature
    glColor4f(0.0,1.0,0.0,1.0);
    for(int i = 0; i < n; i++){
        glBegin(GL_LINE_LOOP);
        glVertex2f(i*l*u-1.0, 0.0);
        glVertex2f(i*l*u-1.0, f( (i*l*u + (i+1)*l*u) / 2 ));
        glVertex2f((i+1)*l*u-1.0, f( (i*l*u + (i+1)*l*u) / 2 ));
        glVertex2f((i+1)*l*u-1.0, 0.0);
        glEnd();
    }
    glColor4f(0.0,1.0,0.0,0.3);
    for(int i = 0; i < n; i++){
        glBegin(GL_POLYGON);
        glVertex2f(i*l*u-1.0, 0.0);
        glVertex2f(i*l*u-1.0, f( (i*l*u + (i+1)*l*u) / 2 ));
        glVertex2f((i+1)*l*u-1.0, f( (i*l*u + (i+1)*l*u) / 2 ));
        glVertex2f((i+1)*l*u-1.0, 0.0);
        glEnd();
    }
    
    // draw composite trapezoid quadrature
    glColor4f(0.0,0.0,0.0,1.0);
    for(int i = 0; i < n; i++){
        glBegin(GL_LINE_LOOP);
        glVertex2f(i*l*u-1.0, 0.0);
        glVertex2f(i*l*u-1.0, f(i*l*u));
        glVertex2f((i+1)*l*u-1.0, f((i+1)*l*u));
        glVertex2f((i+1)*l*u-1.0, 0.0);
        glEnd();
    }
    glColor4f(0.0,0.0,0.0,0.4);
    for(int i = 0; i < n; i++){
        glBegin(GL_POLYGON);
        glVertex2f(i*l*u-1.0, 0.0);
        glVertex2f(i*l*u-1.0, f(i*l*u));
        glVertex2f((i+1)*l*u-1.0, f((i+1)*l*u));
        glVertex2f((i+1)*l*u-1.0, 0.0);
        glEnd();
    }

    glFlush();
    glutSwapBuffers();
}

void
reshape(int width, int height){
    glViewport(0, 0, width, height); 
    glLoadIdentity();
    if(width > height){								
        glScalef(0.2 * (double) height/width, 0.6, 0.4);
    }
    else{												
        glScalef(0.2, 0.6 * (double) width/height, 0.4);
    }
}

void
keyboard(unsigned char key, int x, int y){
    switch(key){
        case 27: exit(EXIT_SUCCESS);
    }
}


/*--------------------------------------------------------------
 * Main function
 * ------------------------------------------------------------*/
int main(int argc, char** argv){
    a = 0.0;
    b = M_PI * 0.5;
    n = 10;
    
    midquad = setup_midpointrule();
    areamid = eval_quadrature(midquad, a, b, (function)f, data);
    areamid_comp = eval_composite_quadrature(midquad, a, b, n, (function)f, data);

    
    trquad = setup_trapezoidalrule();
    areatr = eval_quadrature(trquad, a, b, (function)f, data);
    areatr_comp = eval_composite_quadrature(trquad, a, b, n, (function)f, data);
    
    printf("Integral Midpoint: %f\n",areamid);
    printf("Integral Trapezoidal: %f\n\n",areatr);
    printf("Composite Quadratures for n = %d:\n",n);
    printf("Composite Midpoint: %f\n",areamid_comp);
    printf("Composite Trapezoidal: %f\n",areatr_comp);
    
    
    
    glutInit(&argc, argv);
    glutCreateWindow("Quadrature");
    glutPositionWindow(100, 100);
    glutReshapeWindow(600, 600); 
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutCreateWindow("Composition");
    glutPositionWindow(700, 100);
    glutReshapeWindow(600, 600);
    glutDisplayFunc(display2);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    
//     glutKeyboardFunc(key2); 
//     del_quadrature(midquad);
//     del_quadrature(trquad);
//     del_pquadrature(trquad);
//     glutTimerFunc(50, timer, 0);
    glutMainLoop();

    return EXIT_SUCCESS;
}