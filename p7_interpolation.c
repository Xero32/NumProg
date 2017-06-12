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
int tmr = 0.5;
// int n = 50;
int n;
int fctflag;
int interflag = 0;
// const double lambda = 1.0;
double data[2] = {0.75,1.0};
double p0, p1, p2;
double a,b;
pinterpolation inter0, inter1, inter2;

/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
double
f(double x, double *lambda){
    assert(*lambda);
    return *lambda * exp(-0.5 * pow(x,2));
}

double
g(double x){
    assert(1 + pow(x,2)); 
    return 1 / (1 + pow(x,2));
}

void
setup_transition_points(pinterpolation inter0, pinterpolation inter1, pinterpolation inter2, double a, double b){
    double *x0 = inter0->xi;
    double *x1 = inter1->xi;
    double *x2 = inter2->xi;
    tmr = 0.5;
    for(int i = 0; i < m; i++){
//         x2[i] = (x1[i] - x0[m-i]) * 0.5 + x0[i];
        printf("t: %d\n",tmr);
        printf("x0[%d]: %f, x1[%d]: %f, x2[%d]: %f\n", i,x0[m-i],i,x1[i],i,x2[i]);
    }

//     t += 0.5;
}

/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/
void
display(){
    glutSetWindow(1);
    double *y1 = inter1->f;
    double *x1 = inter1->xi;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glPushMatrix();
    //draw interpolation
    
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.8,0.0);
    for(int i = 0; i < m; i++){
        glVertex2f(x1[i], y1[i] - 0.5);
//         printf("i: %d, xi[i]: %f, y[i]: %f\n",i,x1[i],y1[i]);
    }
    glEnd();

    // draw red squares at interpolation points
    for(int i = 0; i < m; i++){
        glBegin(GL_LINE_LOOP);
        glColor3f(1.0,0.0,0.0);
        glVertex2f(x1[i] - 0.01, y1[i] - 0.51);
        glVertex2f(x1[i] + 0.01, y1[i] - 0.51);
        glVertex2f(x1[i] + 0.01, y1[i] - 0.49);
        glVertex2f(x1[i] - 0.01, y1[i] - 0.49);
        glEnd();    
    }
    
    // draw interpolated polynomial
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    newton_divided_differences(inter1);
    for(int j = 0; j < 300; j++){
        p1 = eval_interpolation_polynomial(inter1, (double)  (j/100.0)-1.0) - 0.5;
        glVertex2f((double)(j/100.0)-1.0,p1);
    }
    glEnd();
    
    double *x2 = inter2->xi;
    double *y2 = inter2->f;
//     setup_transition_points(inter0, inter1, inter2, a,b);
//     eval_interpolated_values(inter2, (function)f, data);
//     newton_divided_differences(inter2);
    // connect interpolation points
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.8,0.9);
    for(int i = 0; i < m; i++){
        glVertex2f(x2[i], y2[i] - 0.5);
    }
    glEnd();
    
    // draw red squares at interpolation points
    for(int i = 0; i < m; i++){
        glBegin(GL_LINE_LOOP);
        glColor3f(1.0,0.0,1.0);
        glVertex2f(x2[i] - 0.01, y2[i] - 0.51);
        glVertex2f(x2[i] + 0.01, y2[i] - 0.51);
        glVertex2f(x2[i] + 0.01, y2[i] - 0.49);
        glVertex2f(x2[i] - 0.01, y2[i] - 0.49);
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
        glScalef((double) height/width, 1.0, 1.0);
    }
    else{												
        glScalef(1.0, (double) width/height, 1.0);
    }
}

void
select_function(int flag){
    printf("flag: %d\n", flag);
    if(!flag){
        eval_interpolated_values(inter1, (function)f, data);
        newton_divided_differences(inter1);
        glutPostRedisplay();
    }else{
        eval_interpolated_values(inter1, (function)g, data);
        newton_divided_differences(inter1);
        glutPostRedisplay(); 
    }
}
/* If necessary use the pictures as orientation,
   but feel free to create something different =)*/

void
key(unsigned char key, int x, int y){
    switch(key){
        case 27:    exit(EXIT_SUCCESS);
        case 's':   if(!interflag){ 
//                         printf("fctflag: %d\n",fctflag);
                        setup_aequidistant_interpolationpoints(inter1, a,b);
                        select_function(!fctflag); 
                        interflag = 1 - interflag;
                        glutPostRedisplay(); break;
                    }else{
//                         printf("fctflag: %d\n",fctflag);
                        setup_chebyshev_interpolationpoints(inter1, a,b);
                        select_function(!fctflag);
                        interflag = 1 - interflag;
                        glutPostRedisplay(); break;
                    }
        case 'f':   if(!fctflag){
//                         printf("fctflag: %d\n",fctflag);
                        select_function(fctflag);
                        fctflag = 1 - fctflag;
                        glutPostRedisplay(); break;
                    }else{
//                         printf("fctflag: %d\n",fctflag);
                        select_function(fctflag);
                        fctflag = 1 -  fctflag;
                        glutPostRedisplay(); break;
                    }
    
         

    }
}

// void
// timer(int val){
//     double *x0 = inter0->xi;
//     double *x1 = inter1->xi;
//     double *x2 = inter2->xi;
// //     for(int i = 0; i < m; i++){
// //         x2[i] = (x1[i] - x0[i]) * t + x0[i];
// //     }
// //     t += 0.05;
//     eval_interpolated_values(inter2, (function)f, data);
//     glutPostRedisplay();
//     glutTimerFunc(20, timer, 0);
// }

void 
transition(){
    glutSetWindow(2);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.0,1.0);
    glPushMatrix();
//     double *x0 = inter0->xi;
//     double *x1 = inter1->xi;
    double *x2 = inter2->xi;
    double *y2 = inter2->f;
    // draw red squares at interpolation points
    for(int i = 0; i < m; i++){
        
        glBegin(GL_LINE_LOOP);
        glColor3f(1.0,0.0,1.0);
        glVertex2f(x2[i] - 0.01, y2[i] - 0.51);
        glVertex2f(x2[i] + 0.01, y2[i] - 0.51);
        glVertex2f(x2[i] + 0.01, y2[i] - 0.49);
        glVertex2f(x2[i] - 0.01, y2[i] - 0.49);
        glEnd();    
    }
    
    
/*    newton_divided_differences(inter1);
    
    glBegin(GL_LINE_STRIP);
    glColor3f(0.9,0.2,0.2);
    for(int i = 0; i < 300; i++){
        x1[i] = (x1[i] - x0[i]) * t + x0[i];
        p1 = eval_interpolation_polynomial(inter1, (double)  (i/100.0)-1.0) - 0.5;
        glVertex2f((double)(i/100.0)-1.0,p1); 
    }
    glEnd();
    glFlush();
    glutSwapBuffers();  */ 
//     t += 0.01;
//     glutTimerFunc(15, transition, 0);  
//     glutPostRedisplay();
}


int
main(int argc, char **argv){
//     assert(n && m);

    a = -0.5;
    b = 0.5;
    inter0 = new_interpolation(m);
    setup_aequidistant_interpolationpoints(inter0, a,b);
    eval_interpolated_values(inter0, (function)g, data);
    newton_divided_differences(inter0);
    
    inter1 = new_interpolation(m);
    setup_chebyshev_interpolationpoints(inter1, a,b);
    eval_interpolated_values(inter1, (function)g, data);
    newton_divided_differences(inter1);

    
    inter2 = new_interpolation(m);
//     setup_transition_points(inter0, inter1, inter2, a,b);
//     eval_interpolated_values(inter2, (function)f, data);
//     newton_divided_differences(inter2);
    
    
    glutInit(&argc, argv);
    glutCreateWindow("Interpolation");
    glutPositionWindow(100, 100);
    glutReshapeWindow(600, 600); 
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    
//     glutCreateWindow("Transition");
//     glutPositionWindow(700, 700);
//     glutReshapeWindow(600, 600);
//     glutReshapeFunc(reshape);
//     glutDisplayFunc(transition);
//     
    glutKeyboardFunc(key);
//     glutKeyboardFunc
//      glutTimerFunc(50, timer, 0);
    glutMainLoop();
    
    
    
    
    
    

    return EXIT_SUCCESS;
}