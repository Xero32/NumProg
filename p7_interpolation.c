/*---------------------------------------------------------------*/
/*		     Numerische Programmierung  	 	 */
/* 	Serie 7 - Approximation von Funktionen (Interpolation) 	 */
/* ------------------------------------------------------------- */
/*	Autoren: Christina Boerst				 */
/*	Versionsnummer:	1					 */
/*---------------------------------------------------------------*/
/*      Bearbeitet: Marko Hollm, Marvin Becker                   */
/*      Version:    1                                            */
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

#define CHKPOS(f) if((err=f) < 0) goto Error
/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
int m = 6;
double t = 0.0;
// int n = 50;
int n;
int fctflag, fctflag2;
int interflag = 0;
double data[1] = {2.9};
double p0, p1;
double a,b;
int err;
pinterpolation inter0, inter1;
int FLAG = 0;
/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
static void
Printhelp(){
    printf("Help:\nThis program displays the interpolation of functions f and g given in the task.\n\n");
    printf("For function f you can adjust the coefficient lambda when starting the program.\n");
    printf("Usage example for Linux:\t'./p7_interpolation 0.5' sets lambda to 0.5\n\n");
    printf("In the 'Interpolation' window you see the function with \nits static interpolation points.\n");
    printf("Here you can switch between functions f and g by pressing 'f'.\nChange interpolation mode from aequidistant to Chebyshev\n");
    printf("or vice versa with 's'.\n\nIn the 'Transition' window you can watch the interpolation nodes convert \nfrom aequidistant to Chebyshev.\n");
    printf("Reset with 'r', switch between functions f and g with 'f'\n\n");
    printf("In both windows you can increase/decrease the number of interpolation nodes \nwith '+' or '-' respectively.\n\n");
    printf("Set the interval boundaries with 'a' and 'b'. \nAfter pressing the respective button, \nyou will have to switch to the shell to put in the value\n");    
}


static void
delete_all(){
    del_interpolation(inter0);
    del_interpolation(inter1);
}

static double
f(double x, double *lambda){
    assert(*lambda);
    return *lambda * exp(-0.5 * pow(x,2));
}

static double
g(double x){
    assert(1 + pow(x,2)); 
    return 1 / (1 + pow(x,2));
}

// to help visualize the transition between interpolation types
// non-instantaneous mapping of one point onto the other
static void
setup_transition_points(pinterpolation inter0, pinterpolation inter1, double step){
    double *x0 = inter0->xi;
    double *x1 = inter1->xi;
    for(int i = 0; i < m+1; i++){
        x0[i] = (x1[i] - x0[i]) * step + x0[i];
    }
}

// select function f or g
static void
select_function(pinterpolation intr, int flag){
    if(!flag){
        eval_interpolated_values(intr, (function)f, data);
        newton_divided_differences(intr);
    }else{
        eval_interpolated_values(intr, (function)g, data);
        newton_divided_differences(intr);
    }
}

// sets up interpolations for further use
static void
setup(int flag){
    inter0 = new_interpolation(m);
    setup_aequidistant_interpolationpoints(inter0, a,b);
    select_function(inter0, !flag);
    
    inter1 = new_interpolation(m);
    setup_chebyshev_interpolationpoints(inter1, a,b);
    select_function(inter1, !flag);
}
/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/
void
timer(int val){
    if(t <= 1.0 && t >= 0.0){
        t += 0.005;
        glutPostRedisplay();
        glutTimerFunc(50, timer, 0);
    }
}

//display function for 'Interpolation' window
void
display(){
    glutSetWindow(1);
    double *y1 = inter1->f;
    double *x1 = inter1->xi;
    double *d1 = inter1->d;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glPushMatrix();
    
    //draw interpolation
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.8,0.0);
    for(int i = 0; i < m+1; i++){
        glVertex2f(x1[i], y1[i]);
        if(FLAG == 1) printf("d[%d] = %f\n",i,d1[i]);
    }
    glEnd();

    // draw red squares at interpolation points
    for(int i = 0; i < m+1; i++){
        glBegin(GL_LINE_LOOP);
        glColor3f(1.0,0.0,0.0);
        glVertex2f(x1[i] - 0.04, y1[i] - 0.04);
        glVertex2f(x1[i] + 0.04, y1[i] - 0.04);
        glVertex2f(x1[i] + 0.04, y1[i] + 0.04);
        glVertex2f(x1[i] - 0.04, y1[i] + 0.04);
        glEnd();    
    }
    
    // draw interpolated polynomial
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    newton_divided_differences(inter1);
    for(int j = 0; j < 500; j++){
        p1 = eval_interpolation_polynomial(inter1, (double)(j/50.0)-5.0);
        glVertex2f((double)(j/50.0)-5.0,p1);
    }
    glEnd();
    
    //draw coordinate system
    glBegin(GL_LINES);
    glColor3f(0.0,0.0,0.0);
    glVertex2f(0.0,-100.0);
    glVertex2f(0.0,100.0);
    glVertex2f(100.0,0.0);
    glVertex2f(-100.0,0.0);
    glVertex2f(1.0,-.1);  
    glVertex2f(1.0,.1);
    glVertex2f(-1.0,-.1);
    glVertex2f(-1.0,.1);
    glVertex2f(.1,1.0);
    glVertex2f(-.1,1.0);
    glVertex2f(.1,-1.0);
    glVertex2f(-.1,-1.0);
    glEnd();
    
    glFlush();
    glutSwapBuffers();
    FLAG++;
}

void
display2(){
    glutSetWindow(2);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    // prepare transitioning nodes
    double *x0 = inter0->xi;
    double *y0 = inter0->f;
    setup_transition_points(inter0,inter1,t);
    select_function(inter0, !fctflag2);
    
     // draw purple crosses at interpolation points
    glLineWidth(2.0);
    for(int i = 0; i < m+1; i++){
        glBegin(GL_LINES);
        glColor3f(1.0,0.0,1.0);
        glVertex2f(x0[i] - 0.05, y0[i] - 0.05);
        glVertex2f(x0[i] + 0.05, y0[i] + 0.05);
        glVertex2f(x0[i] - 0.05, y0[i] + 0.05);
        glVertex2f(x0[i] + 0.05, y0[i] - 0.05);
        glEnd();  
    }
    glLineWidth(1.0);

    // draw interpolated polynomial
    glBegin(GL_LINE_STRIP);
    glColor3f(0.1,0.1,0.1);
    for(int i = 0; i < 500; i++){
        p0 = eval_interpolation_polynomial(inter0, (double) (i/50.0)-5.0);
        glVertex2f((double)(i/50.0)-5.0,p0);
    }
    glEnd();
    
    //draw coordinate system
    glBegin(GL_LINES);
    glColor3f(0.0,0.0,0.0);
    glVertex2f(0.0,-100.0);
    glVertex2f(0.0,100.0);
    glVertex2f(100.0,0.0);
    glVertex2f(-100.0,0.0);
    glVertex2f(1.0,-.1);  
    glVertex2f(1.0,.1);
    glVertex2f(-1.0,-.1);
    glVertex2f(-1.0,.1);
    glVertex2f(.1,1.0);
    glVertex2f(-.1,1.0);
    glVertex2f(.1,-1.0);
    glVertex2f(-.1,-1.0);
    glEnd();
    
    glFlush();
    glutSwapBuffers();
}

void
reshape(int width, int height){
    glViewport(0, 0, width, height); 
    glLoadIdentity();
    if(width > height){								
        glScalef(0.2 * (double) height/width, 0.2, 0.2);
    }
    else{												
        glScalef(0.2, 0.2 * (double) width/height, 0.2);
    }
}

/* If necessary use the pictures as orientation,
   but feel free to create something different =)*/

void
key(unsigned char key, int x, int y){
    switch(key){
        case 27:    exit(EXIT_SUCCESS);
        case 's':   if(!interflag){ 
                        setup_aequidistant_interpolationpoints(inter1, a,b);
                        select_function(inter1, !fctflag); 
                        interflag = 1 - interflag;
                        glutPostRedisplay(); return;
                    }else{
                        setup_chebyshev_interpolationpoints(inter1, a,b);
                        select_function(inter1, !fctflag);
                        interflag = 1 - interflag;
                        glutPostRedisplay(); return;
                    }
        case 'f':   if(!fctflag){
                        select_function(inter1, fctflag);
                        fctflag = 1 - fctflag;
                        glutPostRedisplay(); return;
                    }else{
                        select_function(inter1, fctflag);
                        fctflag = 1 -  fctflag;
                        glutPostRedisplay(); return;
                    }
        case 0x02B: m++; // '+' key;
                    delete_all();
                    setup(fctflag);
                    glutPostRedisplay();
                    return;
        case 0x02D: if(m>1){m--; //'-' key;
                        delete_all();
                        setup(fctflag);
                        glutPostRedisplay();
                        return;
                    }else{
                        printf("You can't decrease the amount of nodes any further. \n");
                        return;	
                    }
        case 'a':   printf("choose left boundary: \n");
                    CHKPOS(scanf("%lf",&a));
                    delete_all();
                    setup(fctflag);
                    glutPostRedisplay();
                    return;
        case 'b':   printf("choose right boundary: \n");
                    CHKPOS(scanf("%lf",&b));
                    delete_all();
                    setup(fctflag);
                    glutPostRedisplay();
                    return;
        case 'h':   Printhelp();
                    glutSetWindow(1);
                    glutIconifyWindow();
                    glutSetWindow(2);
                    glutIconifyWindow();
                    return;
        default:    return;
    }
    Error:
    printf("You have not entered a number.\n");
    return;
}

void
key2(unsigned char key, int x, int y){
    switch(key){
        case 27:    exit(EXIT_SUCCESS);
        case 'r':   t = 0.0;
                    delete_all();
                    setup(fctflag);
                    glutTimerFunc(50,timer,0);
                    return;
        case 'f':   if(!fctflag2){
                        select_function(inter1, fctflag2);
                        fctflag2 = 1 - fctflag2;
                        glutPostRedisplay(); return;
                    }else{
                        select_function(inter1, fctflag2);
                        fctflag2 = 1 -  fctflag2;
                        glutPostRedisplay(); return;
                    }
        case 0x02B: m++; // '+' key;
                    delete_all();
                    setup(fctflag);  
                    glutPostRedisplay();
                    return;
        case 0x02D: if(m>1){m--; //'-' key;
                        delete_all();
                        setup(fctflag);
                        glutPostRedisplay();
                        return;
                    }else{
                        printf("You can't decrease the amount of nodes any further. \n");
                        return;	
                    }	
       case 'a':    printf("choose left boundary: \n");
                    CHKPOS(scanf("%lf",&a));
                    delete_all();
                    setup(fctflag);
                    glutPostRedisplay();
                    return;
        case 'b':   printf("choose right boundary: \n");
                    CHKPOS(scanf("%lf",&b));
                    delete_all();
                    setup(fctflag);
                    glutPostRedisplay();
                    return;
        case 'h':   Printhelp();
                    glutSetWindow(1);
                    glutIconifyWindow();
                    glutSetWindow(2);
                    glutIconifyWindow();
                    return;
        default:    return;
    }
    Error:
    printf("You have not entered a number.\n");
    return;
}

int
main(int argc, char **argv){
//     assert(n && m);
    if(argc > 1){
        data[0] = atof(argv[1]);
        printf("lambda: %f\n", data[0]);
    }
    a = 2.1;
    b = -2.1;
        
    setup(fctflag);   
    
    glutInit(&argc, argv);
    glutCreateWindow("Interpolation");
    glutPositionWindow(100, 100);
    glutReshapeWindow(600, 600); 
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key);
    glutCreateWindow("Transition");
    glutPositionWindow(700, 100);
    glutReshapeWindow(600, 600);
    glutReshapeFunc(reshape);
    glutDisplayFunc(display2);
    
    glutKeyboardFunc(key2); 
    
    glutTimerFunc(50, timer, 0);
    glutMainLoop();

    return EXIT_SUCCESS;
}