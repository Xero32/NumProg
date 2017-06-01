/*---------------------------------------------------------------*/
/*		     Numerische Programmierung 	 	 	 */
/* 	Serie 6 - Visualisierung von Oberflaechentriangulationen     */
/* ------------------------------------------------------------- */
/*	Autoren: 		Ove Hansen, Tim Drevelow			 		 */
/*	Versionsnummer:			1.0									 */
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

							
psurface3d gr; 

int mouseX,mouseY;
int windowWidth,windowHeight;
char inputmode='r';
float angle_x,angle_y;
float tx,ty,tz;
float zoom=1.0;
int wire=0;


/* Translation */

static void
translate(double x, double y, double z){
 
	double translate[16]={1,0,0,0,0,1,0,0,0,0,1,0,x,y,z,1};
	glMultMatrixd(translate);
  
}

/* Rotation around x-axis */
static void
rotate_x(double alpha){
		
	double rotate_x[16]={1,0,0,0,0,cos(alpha/720*M_PI),sin(alpha/720*M_PI),0,0,-sin(alpha/720*M_PI),cos(alpha/720*M_PI),0,0,0,0,1};
	glMultMatrixd(rotate_x);
  
}

/* Rotation around y-axis */
static void
rotate_y(double alpha){
	
	
	double rotate_y[16]={cos(alpha/720*M_PI),0,-sin(alpha/720*M_PI),0,0,1,0,0,sin(alpha/720*M_PI),0,cos(alpha/720*M_PI),0,0,0,0,1};
	glMultMatrixd(rotate_y);
}

static void
display_mesh_full(){
	
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glEnable(GL_DEPTH_TEST);

	
	glEnable(GL_LIGHTING); 
	glEnable(GL_LIGHT0);
	
	float light_position[4] = { 5.0, 5.0, 5.0, 1.0 }; 
	float light_ambient[4] = { 0.1, 0.1, 0.1, 1.0 }; 
	float light_diffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position); 
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient); 
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse); 

	

	int triangles=gr->triangles;
	int (*t)[3]=gr->t;
	real (*x)[3]=gr->x;
	real (*n)[3]=gr->n;
	
	glPushMatrix();
	
	translate(tx,ty,tz);
	rotate_x(angle_x);
	rotate_y(angle_y);
	
	glBegin(GL_TRIANGLES);
	
	float material_ambient[4] = { 1.0, 0.0, 0.0, 0.0 }; 
	float material_diffuse[4] = { 1.0, 0.0, 0.0, 0.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient); 
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);

	for(int i=0; i<triangles; i++){
		glNormal3f(n[i][0],n[i][1],n[i][2]);
		for(int j=0; j<3; j++){
			glVertex3f(x[t[i][j]][0],x[t[i][j]][1],x[t[i][j]][2]);
		}
	}
	glEnd();
	glPopMatrix();	
}

static void
display_mesh_wire(){
	
	glEnable(GL_DEPTH_TEST);
	
	glDisable(GL_LIGHTING); 
	glDisable(GL_LIGHT0);

	int edges=gr->edges;
	int (*e)[2]=gr->e;
	real (*x)[3]=gr->x;
	
	glPushMatrix();
	
	translate(tx,ty,tz);
	rotate_x(angle_x);
	rotate_y(angle_y);
	
	glBegin(GL_LINES);

	for(int i=0; i<edges; i++){
		for(int j=0; j<2; j++){
			glVertex3f(x[e[i][j]][0],x[e[i][j]][1],x[e[i][j]][2]);
		}
	}
	glEnd();
	glPopMatrix();		
}
  
 /* Drawing complete surface triangulation */
static void
display_mesh(){
	
	glClearColor(0.3f, 0.3f, 0.3f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glPushMatrix();
	glLoadIdentity();
	gluPerspective(35.0, 1.0, 1.0*zoom, 30.0*zoom);
	glTranslatef(0.0, 0.0, -10.0);
	glPopMatrix();
	
	if(wire==0)display_mesh_full();
	if(wire==1)display_mesh_wire();
	
	glFlush();
    glutSwapBuffers();
}
  
  /* Reshape of a triangulation */
static void
reshape_mesh(int width, int height){
	
	windowWidth=width;
	windowHeight=height;
    glViewport(0, 0, width, height); 
	glLoadIdentity();
	
	if(width > height){								
       glScalef((double) height/width, 1.0, 1.0);
    }
    else{												
       glScalef(1.0, (double) width/height, 1.0);
    }
    gluPerspective(35.0, 1.0, 1.0, 30.0);
 	glTranslatef(0.0, 0.0, -10.0);

}

  /* Mouse movement on a triangulation */
static void
mouse_mesh(int button, int state, int position_y, int position_x){
	mouseX=position_x;
	mouseY=position_y;
}

  /* Motion on a triangulation */
static void
motion_mesh(int position_y, int position_x){
	if(inputmode=='r'){
		angle_y+=(float)(position_y-mouseY)/windowHeight*360;
		angle_x+=(float)(position_x-mouseX)/windowWidth*360;
	} else if(inputmode=='t'){
		ty-=(float)(position_x-mouseX)/windowWidth*10;
		tx+=(float)(position_y-mouseY)/windowHeight*10;
	}
	mouseY=position_y;
	mouseX=position_x;
	glutPostRedisplay();
}
 
  /* Key input for triangulation */
static void
key_mesh(unsigned char key, int x, int y){
	
	if(key=='r')inputmode=key;
	
	if(key=='t')inputmode=key;
	
	if(key==' '){
		angle_x=0.0;
		angle_y=0.0;
		tx=0.0;
		ty=0.0;
		tz=0.0;
		inputmode='r';
		glutPostRedisplay();
	}
	if(key=='w')wire = 1 - wire;
	if(key=='q')zoom *= 0.9;
	if(key=='a')zoom /= 0.9;
  
}
   
  
/* last but not least the main function */   
int
main(int argc,char **argv){
  
	psurface3d sur;
	//int i;  I really don't know what this is for TODO
  
	/* Reading mesh */
    if(argc > 1){
      sur = read_surface3d(argv[1]); 
  	}
    else{
      printf("No input file!\n");
      return EXIT_SUCCESS;
    }
	
	prepare_surface3d(sur);
	gr=sur;
  
    /*-------------------------------------------------
	*	 Initalize GUI
	*--------------------------------------------------*/
    glutInit(&argc, argv);
	glutCreateWindow("Mesh");
	
	glutPositionWindow(150, 100);
	glutReshapeWindow(1000, 700);
	
	/*callback functions*/
	glutReshapeFunc(reshape_mesh);
	glutDisplayFunc(display_mesh);
    glutKeyboardFunc(key_mesh);
    glutMotionFunc(motion_mesh);
	glutMouseFunc(mouse_mesh);
	
	glutMainLoop();
   return EXIT_SUCCESS;
}
