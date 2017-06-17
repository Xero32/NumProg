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
    #define MAXACC 100
    #define CHKPOS(f) if((err=f) < 0) goto Error
/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
//TODO: composite trapezoidal gives wrong result for n = 1!
double a, b, areamid, areatr;
double areamid_comp, areatr_comp;
int m, n;
int quadflag;
pquadrature midquad, trquad;
double data[1] = {0.0};
int err;
char _areamid[50] = {'\0'}; char _areatr[50] = {'\0'}; char _n[50] = {'\0'};
char _areamid_comp[50] = {'\0'}; char _areatr_comp[50] = {'\0'};

static void 
Printhelp(){
    printf("This program shows the single midpoint (green area) and \ntrapezoidal (grey area) quadrature in window 'Quadrature'.\n");
    printf("The corresponding area values are displayed, if the interval \nis chosen to be close to [0,pi/2].\n\n");
    printf("In window 'Composition' you can see the composite quadratures \nof the same kinds, with the order 'n' and areas \non screen for appropriate intervals.\n\n");
    printf("Optionally you can set 'n' as an input at the start of the program. Afterwards\n");
    printf("increase or decrease 'n' with the '+' or '-' key and watch \nthe calculated area converge.\n");
    printf("The areas are given in the command prompt as well.\n");
    
}
/* ------------------------------------------------------------
* Example function
*-------------------------------------------------------------*/
static void
swap(double *a, double *b){
        double c = *a;
        *a = *b;
        *b = c;    
}

static double
f(double x){
    return 5.0 * exp(2.0 * x) * cos(x) / (exp(M_PI)-2.0);
}

double
fct_max(function f, void *data, double a, double b){
    if(a > b) swap(&a, &b);
    double max = 0.0;
    double l = (b-a) * 0.5;
    double u = 1.0 / (MAXACC);
    
    for(int i = 0; i < l * 10.0 * MAXACC; i++){
        if((double)i*l*u+a > b) return max;
        max = ( fabs((double)f(i*l*u+a,data)) > max) ? fabs(f((double)i*l*u+a,data)) : max;
    }
    return max;
}

static void 
output(float x, float y, float r, float g, float b, void *font, char *string){
  glColor3f( r, g, b );
  glRasterPos2f(x, y);
  int len, i;
  len = 50;
  for (i = 0; i < len; i++) {
    glutBitmapCharacter(font, string[i]);
  }
}

static void
setup(pquadrature quad, int flag, double *area, double *area_comp, function f, void *data){
    switch(flag){
        case 1: quad = setup_midpointrule(); break;
        case 2: quad = setup_trapezoidalrule(); break;
    }
    if(a > b) swap(&a, &b);
    *area = eval_quadrature(quad, a, b, f, data);
    if(n > 1) *area_comp = eval_composite_quadrature(quad, a, b, n, f, data);
    else *area_comp = *area; //eval_composite_quadrature(quad, a, b, 1, f, data);  //TODO Problem!
    //     else *area_comp = eval_quadrature(quad, a, b, f, data);
    del_quadrature(quad);
}

/* ------------------------------------------------------------
* GLUT functions
*-------------------------------------------------------------*/
static void
draw_coord(){
    glColor3f(0.0,0.0,0.0);
    glVertex2f(0.0,-1000.0);    //y-axis
    glVertex2f(0.0,1000.0);     //y-axis
    glVertex2f(1000.0,0.0);      //x-axis
    glVertex2f(-1000.0,0.0);     //x-axis
    glVertex2f(-0.5,-.03);      //x-mark
    glVertex2f(-0.5,.03);
    glVertex2f(0.0,-.07);       
    glVertex2f(0.0,.07);
    glVertex2f(0.5,-.03);
    glVertex2f(0.5,.03);
    glVertex2f(1.0,-.07);
    glVertex2f(1.0,.07);
    glVertex2f(1.5,-.03);
    glVertex2f(1.5,.03);
    glVertex2f(5.0,-100.0);
    glVertex2f(5.0,100.0);

    glVertex2f(0.1,1.0);       //y-mark
    glVertex2f(-0.1,1.0);  
    glVertex2f(-0.05,0.5);
    glVertex2f(0.05,0.5);
    glVertex2f(-0.05,-0.5);
    glVertex2f(0.05,-0.5);
    glVertex2f(0.1,-1.0);      //y-mark
    glVertex2f(-0.1,-1.0);
}

static void
draw_marks(double loc, double wid){
    glVertex2f(loc, wid*0.5);
    glVertex2f(loc,-wid*0.5);  
}

static void
display(){
    if(a > b) swap(&a, &b);
    glutSetWindow(1);
    double aa = a;
    double bb = b;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glEnable (GL_BLEND); 
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    
    // draw function
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    for(int j = 0; j < 500; j++){
        glVertex2f((double)(j/50.0)+a-1.0,f((double)(j/50.0)+a-1));
    }
    glEnd();
    
    // draw single midpoint quadrature area
    glColor4f(0.0,1.0,0.0,0.4);
    glBegin(GL_POLYGON);
    glVertex2f(aa,0.0);
    glVertex2f(aa,f((a+b)*0.5));
    glVertex2f(bb,f((a+b)*0.5));
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
    draw_marks(aa,.2);
    draw_marks(bb,.2);
    glEnd();
    glLineWidth(1.0);
    
    // Omit display of area for vastly bigger or smaller boundaries, 
    // because scaling with max(f) ruins arrangement of text.
    // Show information on screen only for this particular function 
    // near the given interval [0,pi/2].
    // It is not meant for generalization, just a little treat.
    if(a > -1.0 && a < 1.0 && b > M_PI * 0.5 - 0.2 && b <= 1.8){
        snprintf(_areamid, 50, "%f", areamid);
        snprintf(_areatr, 50, "%f", areatr);
        
        output(1.0, 0.5, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, _areamid);
        output(0.2, 0.05, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, _areatr);
    }
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}

static void
display2(){
    if(a > b) swap(&a, &b);
    glutSetWindow(2);
    double aa = a;
    double bb = b;
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.9,0.9,0.9,1.0);
    glEnable (GL_BLEND); 
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPushMatrix();
    
     // draw function
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,0.9);
    for(int j = 0; j < 500; j++){
        glVertex2f((double)(j/50.0)+a-1.0,f((double)(j/50.0)+a-1));
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
    draw_marks(aa,.2);
    draw_marks(bb,.2);
    glEnd();
    glLineWidth(1.0);
    
    double u = 1.0 / (double)n;
    double l = b-a;
    
    //draw composite midpoint quadrature
    glColor4f(0.0,1.0,0.0,1.0);
    for(int i = 0; i < n; i++){
        glBegin(GL_LINE_LOOP);
        glVertex2f(i*l*u+a, 0.0);
        glVertex2f(i*l*u+a, f( (i*l*u + (i+1)*l*u + 2.0 * a) * 0.5 ));
        glVertex2f((i+1)*l*u+a, f( (i*l*u + (i+1)*l*u + 2.0 * a) * 0.5 ));
        glVertex2f((i+1)*l*u+a, 0.0);
        glEnd();
    }
    glColor4f(0.0,1.0,0.0,0.3);
    for(int i = 0; i < n; i++){
        glBegin(GL_POLYGON);
        glVertex2f(i*l*u+a, 0.0);
        glVertex2f(i*l*u+a, f( (i*l*u + (i+1)*l*u + 2.0 * a) * 0.5 ));
        glVertex2f((i+1)*l*u+a, f( (i*l*u + (i+1)*l*u + 2.0 * a) * 0.5 ));
        glVertex2f((i+1)*l*u+a, 0.0);
        glEnd();
    }
    
    // draw composite trapezoid quadrature
    glColor4f(0.0,0.0,0.0,1.0);
    for(int i = 0; i < n; i++){
        glBegin(GL_LINE_LOOP);
        glVertex2f(i*l*u+a, 0.0);
        glVertex2f(i*l*u+a, f(i*l*u+a));
        glVertex2f((i+1)*l*u+a, f((i+1)*l*u+a));
        glVertex2f((i+1)*l*u+a, 0.0);
        glEnd();
    }
    glColor4f(0.0,0.0,0.0,0.4);
    for(int i = 0; i < n; i++){
        glBegin(GL_POLYGON);
        glVertex2f(i*l*u+a, 0.0);
        glVertex2f(i*l*u+a, f(i*l*u+a));
        glVertex2f((i+1)*l*u+a, f((i+1)*l*u+a));
        glVertex2f((i+1)*l*u+a, 0.0);
        glEnd();
    }

    // Omit display of area for vastly bigger or smaller boundaries, 
    // because scaling with max(f) ruins arrangement of text.
    // Show information on screen only for this particular function 
    // near the given interval [0,pi/2].
    // It is not meant for generalization, just a little treat.
    if(a > -1.0 && a < 1.0 && b > M_PI * 0.5 - 0.2 && b <= 1.8){
        snprintf(_areamid_comp, 50, "%f", areamid_comp);
        char txt1[50] = "Composite Midpoint Area:";
        output(a+0.02, -0.4, 0.0,0.0,0.0,GLUT_BITMAP_HELVETICA_12, txt1);
        output(a+0.05, -0.5, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, _areamid_comp);
        glColor4f(0.0,1.0,0.0,0.3);
        glBegin(GL_POLYGON);
        glVertex2f(a+0.02,-0.55);
        glVertex2f(a+0.02,-0.45);
        glVertex2f(a+0.38,-0.45);
        glVertex2f(a+0.38,-0.55);
        glEnd();
        
        snprintf(_areatr_comp, 50, "%f", areatr_comp);
        char txt2[50] = "Composite Trapezoid Area:";
        output(a+0.02, -0.1, 0.0, 0.0, 0.0,GLUT_BITMAP_HELVETICA_12, txt2);
        output(a+0.05, -0.2, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, _areatr_comp);
        glColor4f(0.0,0.0,0.0,0.3);
        glBegin(GL_POLYGON);
        glVertex2f(a+0.02,-0.25);
        glVertex2f(a+0.02,-0.15);
        glVertex2f(a+0.38,-0.15);
        glVertex2f(a+0.38,-0.25);
        glEnd();
        
        snprintf(_n, 50, "%d", n);    
        char txtn[50] = "n = ";
        output(0.25, 1.0, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, txtn);
        output(0.34, 1.0, 0.0, 0.0, 0.0, GLUT_BITMAP_HELVETICA_12, _n);
    }
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}

static void
reshape(int width, int height){
    if(a > b) swap(&a, &b);
    glViewport(0, 0, width, height); 
    glLoadIdentity();
    
    if(width > height){								
        glScalef(0.8 / (b+0.2-a) * 2.0 * (double) height/width, 1.0 / fct_max((function)f,data,a,b) * 0.8, 1.0);
    }
    else{												
        glScalef(0.8 / (b+0.2-a) * 2.0, 1.0 / fct_max((function)f,data,a,b) * 0.8 * (double) width/height, 1.0);
    }
    glTranslatef(-(b+a)/2.0,0.0,0.0);
}

static void
keyboard(unsigned char key, int x, int y){
    switch(key){
        case 27: 
            exit(EXIT_SUCCESS);
        case 0x02B: // '+' key
            glutSetWindow(2); 
            n++; 
            setup(midquad, 1, &areamid, &areamid_comp, (function)f, data);
            setup(trquad, 2, &areatr, &areatr_comp, (function)f, data);
            printf("Composite Quadratures for n = %d:\n",n);
            printf("Composite Midpoint Area: %f\n",areamid_comp);
            printf("Composite Trapezoidal Area: %f\n\n",areatr_comp);
            glutPostRedisplay(); 
            return;
        case 0x02D: // '-' key
            if(n > 1){
                glutSetWindow(2);
                n--;
                setup(midquad, 1, &areamid, &areamid_comp, (function)f, data);
                setup(trquad, 2, &areatr, &areatr_comp, (function)f, data);
                printf("Composite Quadratures for n = %d:\n",n);
                printf("Composite Midpoint Area: %f\n",areamid_comp);
                printf("Composite Trapezoidal Area: %f\n\n",areatr_comp);
                glutPostRedisplay(); 
                return;
            }else{
                printf("You cannot lower n any further.\n");
                return;
            }
//         case 'a':   printf("choose left boundary: \n");
//                     CHKPOS(scanf("%lf",&a));
//                     
//                     glutPostRedisplay();
//                     return;
//         case 'b':   printf("choose right boundary: \n");
//                     CHKPOS(scanf("%lf",&b));
//                     glutPostRedisplay();
//                     return;
        case 'h':   Printhelp();
                    glutSetWindow(1);
                    glutIconifyWindow();
                    glutSetWindow(2);
                    glutIconifyWindow();
                    return;
        default:    return;
    }
//     Error:
//     printf("You have not entered a number.\n");
    return;
}


/*--------------------------------------------------------------
 * Main function
 * ------------------------------------------------------------*/
int 
main(int argc, char** argv){
    a = 0.0;
    b = M_PI * 0.5;
    n = 11;
    
    if(argc > 1) n = atof(argv[1]);
    setup(midquad, 1, &areamid, &areamid_comp, (function)f, data);
    setup(trquad, 2, &areatr, &areatr_comp, (function)f, data);

    printf("Composite Quadratures for n = %d:\n",n);
    printf("Composite Midpoint Area: %f\n",areamid_comp);
    printf("Composite TrapezoidalArea: %f\n\n",areatr_comp);
    printf("Midpoint Integral Area: %f\n",areamid);
    printf("Trapezoidal Integral Area: %f\n\n",areatr);
    
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
    glutMainLoop();
    
    return EXIT_SUCCESS;
}