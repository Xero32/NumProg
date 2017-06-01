/*---------------------------------------------------------------*/
/*		     Numerische Programmierung 	 	 	 */
/* 	Serie 6 - Visualisierung von Oberflaechentriangulationen */
/* ------------------------------------------------------------- */
/*	Autoren: 				 		 */
/*	Versionsnummer:						 */
/*---------------------------------------------------------------*/

		
#include <stdio.h>		
#include <stdlib.h> 
#include <math.h>
#include <GL/freeglut.h>
#include <GL/glut.h>
  
#include "surface3d.h"	

#ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif
  
 

/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
							
psurface3d sur0,sur1; 
real x0 = -2.0, y0 = 0.0, z0 = -3.0; //coordinates of mesh
real x1 = 2.0, y1 = 2.0, z1 = -10.0; //coordinates of solid
real x = 0.0, y = 0.0, z = 0.0;
// real xtr,ytr,ztr; // x-translation, y-translation, z-translation aka zoom
real rotatex,rotatey;
int k = 0;
/* Translation */

static void
translate(double x, double y, double z){
    glLoadIdentity();
    glPushMatrix();
    
//     glTranslatef(x,y,z);
    
    double tr[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1};
    glMultMatrixd(tr);
    glutPostRedisplay();
}

/* Rotation around x-axis */
// static void
// rotate_x(double alpha){
// 
//   /* ---------------------------------------------- */ 
//   /*                                                */
//   /* T T T T T     O O       D D           O O      */
//   /*     T        O   O      D   D        O   O     */
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T        O   O      D   D        O   O     */
//   /*     T         O O       D D           O O      */
//   /*                                                */ 
//   /* ---------------------------------------------- */
//   
// }

/* Rotation around y-axis */
// static void
// rotate_y(double alpha){
// 
//   /* ---------------------------------------------- */ 
//   /*                                                */
//   /* T T T T T     O O       D D           O O      */
//   /*     T        O   O      D   D        O   O     */
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T        O   O      D   D        O   O     */
//   /*     T         O O       D D           O O      */
//   /*                                                */ 
//   /* ---------------------------------------------- */ 	
//   
// }

  
 /* Drawing complete surface triangulation */
static void
display_mesh(){
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    
    int edges = sur0->edges;
    int (*e)[2] = sur0->e;
    real (*x)[3] = sur0->x;
    
    glClearColor(0.0,0.0,0.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
    
    
//     glPopMatrix();
    glBegin(GL_LINES);
    glMatrixMode(GL_PROJECTION);
//     glPushMatrix();
    glTranslatef(x0,y0,z0);
    glColor3f(1.0,1.0,1.0);
    
    glPushMatrix();
    printf("counter: %d\n",k++);
//     glTranslatef(-2,1.5,0.5);
//     translate(xtr,ytr,ztr);
    for(int i = 0; i < edges; i++){
        for(int j = 0; j < 2; j++){
            glVertex3f(x[e[i][j]][0], x[e[i][j]][1], x[e[i][j]][2]);
        }
    }
    glLoadIdentity();
    glEnd();
    
    glPushMatrix();
    glTranslatef(x1,y1,z1);
//     glTranslatef(2.0,0.0,-5.0);
    glColor3f(0.7,0.2,0.0);
    
    int triangles = sur1->triangles;
    int (*t)[3] = sur1->t;
    real (*v)[3] = sur1->x;
    real(*n)[3] = sur1->n;
    glPushMatrix();
    
    glBegin(GL_TRIANGLES);
    for(int i = 0; i < triangles; i++){
        glNormal3f(n[i][0],n[i][1],n[i][2]);
        for(int j = 0; j < 3; j++){
            glVertex3f(v[t[i][j]][0], v[t[i][j]][1], v[t[i][j]][2]);
        }
    }
    glLoadIdentity();
    glEnd();
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}
  
  /* Reshape of a triangulation */
static void
reshape_mesh(int width, int height){
GLfloat p[4],a[4],d[4];

    glViewport(0, 0, width, height); 
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
        /*
        if(width > height){								
        glScalef((double) height/width, 1.0, 1.0);
        }
        else{												
        glScalef(1.0, (double) width/height, 1.0);
        }

        */
    if(width > height)
        glViewport((width-height)/2, 0, height, height);
    else
        glViewport(0, (height-width)/2, width, width);
    
    printf("RESHAPE\n");
    gluPerspective(35.0, 1.0, 1.0, 30.0);
    glTranslatef(x=0.0, y=0.0, z=-5.0);
//     translate(xtr,ytr,ztr);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();	

    glEnable(GL_LIGHTING);	

    glEnable(GL_LIGHT0);

    p[0] = -2.0; p[1] = 2.0; p[2] = 2.0; p[3] = 1.0;
    glLightfv(GL_LIGHT0, GL_POSITION, p); 	

    a[0] = 0.6; a[1] = 0.6; a[2] = 0.6; a[3] = 1.0;	
    glLightfv(GL_LIGHT0, GL_AMBIENT,  a);

    d[0] = 1.0; d[1] = 1.0; d[2] = 1.0; d[3] = 1.0;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, d); 

    glEnable(GL_COLOR_MATERIAL);   					
    glMaterialfv(GL_FRONT, GL_AMBIENT, a);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, d);

    glEnable(GL_DEPTH_TEST);
}

  /* Mouse movement on a triangulation */

//   static void
// mouse_mesh(int button, int state, int position_y, int position_x){
// 
//   /* ---------------------------------------------- */ 
//   /*                                                */
//   /* T T T T T     O O       D D           O O      */
//   /*     T        O   O      D   D        O   O     */
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T       O     O     D     D     O     O    */ 
//   /*     T        O   O      D   D        O   O     */
//   /*     T         O O       D D           O O      */
//   /*                                                */ 
//   /* ---------------------------------------------- */ 	
//  
// }

  /* Motion on a triangulation */
// static void
// motion_mesh(int position_y, int position_x){
//  
// }
 
  /* Key input for triangulation */
static void
key_mesh(unsigned char key, int x, int y){

    switch (key){
        case 'd': translate(x0+=0.05,y0,z0); break;
        case 'a': translate(x0-=0.05,y0,z0); break;
        case 'w': translate(x0,y0+=0.05,z0); break;
        case 's': translate(x0,y0-=0.05,z0); break;
        case 0x02B: translate(x0,y0,z0+=0.05); printf("ZOOM IN\n");break;
        case 0x02D: translate(x0,y0,z0-=0.05); printf("ZOOM OUT\n");break;
        
//         case 'd': x0+=0.1; printf("xtr: %f\n",x0); glutSwapBuffers(); glutPostRedisplay();break;
//         case 'a': x0-=0.1; printf("xtr: %f\n",x0); glutSwapBuffers(); glutPostRedisplay();break;
//         case 's': y0-=0.1; printf("ytr: %f\n",y0); glutSwapBuffers(); glutPostRedisplay();break;
//         case 'w': y0+=0.1; printf("ytr: %f\n",y0); glutSwapBuffers(); glutPostRedisplay();break;
//         case 0x2B: z0+=0.1; printf("ztr: %f\n",z0);glutSwapBuffers(); glutPostRedisplay();break; // '+' key
//         case 0x2D: z0-=0.1; printf("ztr: %f\n",z0);glutSwapBuffers(); glutPostRedisplay();break; // '-' key
        case 27: exit(EXIT_SUCCESS);
    }
}

  
/* last but not least the main function */   
int
main(int argc,char **argv){
  
//   psurface3d sur;
//   int i;
  
  /* Reading mesh */
  if(argc > 1){
    sur0 = read_surface3d(argv[1]); 
    sur1 = read_surface3d(argv[2]);
	}
  else{
    printf("No input file!\n");
    return EXIT_SUCCESS;
  }
 
	glutInit(&argc, argv);
	glutCreateWindow("Triangulation");
	glutPositionWindow(450, 400);
	glutReshapeWindow(900, 600);
	
	glutReshapeFunc(reshape_mesh);
	glutDisplayFunc(display_mesh);

// 	glutMouseFunc(mouse);

//     glutMotionFunc(motion_mesh);

    glutKeyboardFunc(key_mesh);

//     glutTimerFunc(50, timer, 0);
	
	glutMainLoop();
  return EXIT_SUCCESS;
}
