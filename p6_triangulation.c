/*---------------------------------------------------------------*/
/*		     Numerische Programmierung 	 	 	 */
/* 	Serie 6 - Visualisierung von Oberflaechentriangulationen */
/* ------------------------------------------------------------- */
/*	Autoren: 	Marko Hollm, Marvin Becker 		 */
/*	Versionsnummer:  1                                       */
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
		
#define RDSP glutPostRedisplay();
psurface3d sur0;
real x0 = 0.0, y0 = 0.0, z0 = 0.0; // translation coordinates
real xx = 0.0, yy = 1.0, zz = -12.0; // coordinates of body
real rx,ry, anglex, angley;
int currentx, currenty;
int Height, Width;
int rotateflag;
char meshflag;
double zoom = 1.0;
/* Translation */

static void
Printhelp(){
    printf("\nUsage Info:\n");
    printf("After inputfile, you can denote simulation mode:\n\t'm': mesh\t'b': solid\n");
    printf("Switch between solid and mesh simulation mode by pressing 'm'\n\n");
    printf("For easier use, the local coordinate system is shown. \n");
    printf("Here the red axis corresponds to the x-axis, the green one to the y-axis \nand the blue one to the z-axis.\n");
    printf("Use Keyboard to move, zoom, and rotate objects.\n\n");
    printf("Use 'w', 'a', 's', 'd' to move object along its local x- and y-axes.\n");
    printf("with 'x' and 'y' you can translate the object along its local z-axis.\n\nZoom in with '+', zoom out with '-', or with mouse wheel.\n\n");
    printf("Rotate around the x axis by pressing 'r' and 'f',\n");
    printf("Rotate around the y axis by pressing 'q' and 'e'.\n\n");
    printf("You can also rotate the body with your mouse while \npressing and holding any mouse button.\n\n");
    printf("Reset body with 'c'.\n");
    printf("Close with 'esc'.\n\n");
}


static void
translate(double x, double y, double z){
    double tr[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1};
    glMultMatrixd(tr);
}

/* Rotation around x-axis */
static void
rotate_x(double alpha){
    double rotx[16] = {1,0,0,0, 0,cos(alpha),sin(alpha),0, 0,-sin(alpha),cos(alpha),0, 0,0,0,1};
    glMultMatrixd(rotx);    
}

/* Rotation around y-axis */
static void
rotate_y(double alpha){
    double roty[16] = {cos(alpha),0,-sin(alpha),0, 0,1,0,0, sin(alpha),0,cos(alpha),0, 0,0,0,1};
    glMultMatrixd(roty);
}


  
 /* Drawing complete surface triangulation */
static void
display_mesh(){
    // prepare window
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glClearColor(0.9,0.9,0.9,1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
    glPushMatrix();
    // light
    glEnable(GL_LIGHTING); 
    glEnable(GL_LIGHT0);
    float light_position[4] = { 2.0, -3.0, 4.0, 0.5 }; 
    float light_ambient[4] = { 0.7, 0.7, 0.7, 0.8 }; 
    float light_diffuse[4] = { 0.8, 0.8, 0.8, 0.7 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position); 
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient); 
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse); 
    glEnable(GL_LIGHT1);
    float light_position1[4] = {1.0, -1.0, -4.0, 1.0};
    float light_ambient1[4] = { 0.1, 0.1, 0.1, 0.8 }; 
    float light_diffuse1[4] = { 0.5, 0.5, 0.5, 0.7 };
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1); 
    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1); 
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1); 
    
    // scaling
    if(Width > Height){								
        glScalef((double) zoom*Height/Width, zoom, zoom);
    }else{												
        glScalef(zoom, (double) zoom*Width/Height, zoom);
    }
    glPopMatrix();

    // draw body
    if(meshflag == 'm'){ /* mesh */
        glColor3f(0.3,0.2,0.4);
        glLineWidth(1.3);
        int edges = sur0->edges;
        int (*e)[2] = sur0->e;
        real (*x)[3] = sur0->x;
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHT1);
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);

        for(int i = 0; i < edges; i++){
            for(int j = 0; j < 2; j++){
                glVertex3f(x[e[i][j]][0], x[e[i][j]][1], x[e[i][j]][2]);
            }
        }
    }else if(meshflag == 'b'){ /* solid */
        int triangles=sur0->triangles;
	int (*t)[3]=sur0->t;
	real (*x)[3]=sur0->x;
	real (*n)[3]=sur0->n;
        glBegin(GL_TRIANGLES);
	float material_ambient[4] = { 0.3, 0.2, 0.4, 0.0 }; 
	float material_diffuse[4] = { 0.6, 0.4, 0.8, 0.0 };
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient); 
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);

	for(int i=0; i<triangles; i++){
		glNormal3f(n[i][0],n[i][1],n[i][2]);
		for(int j=0; j<3; j++){
			glVertex3f(x[t[i][j]][0],x[t[i][j]][1],x[t[i][j]][2]);
		}
	}
    }
    glLoadIdentity();
    glEnd();
    
    // draw coordinate system
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
    glDisable(GL_LIGHTING);
    glLineWidth(1.5);
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(2.0,0.0,0.0);
    glColor3f(0.0,1.0,0.0);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(0.0,2.0,0.0);
    glColor3f(0.0,0.0,1.0);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(0.0,0.0,2.0);
    glEnd();
    
    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}
  
  /* Reshape of a triangulation */
static void
reshape_mesh(int width, int height){
    Width = width;
    Height = height;
    glViewport(0, 0, width, height); 
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
        
        if(width > height){								
        glScalef((double) height/width, 1.0, 1.0);
        }
        else{												
        glScalef(1.0, (double) width/height, 1.0);
        }

    gluPerspective(35.0, 1.0, 1.0, 30.0);
    glTranslatef(xx,yy,zz);

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glEnable(GL_DEPTH_TEST);
}

  /* Mouse movement on a triangulation */
static void
mouse_mesh(int button, int state, int position_y, int position_x){
    currentx = position_x;
    currenty = position_y;
    // mouse wheel functions
    if(button == 3){
        glScalef(1.01,1.01,1.01); 
        RDSP;
    }else if(button == 4){
        glScalef(0.99,0.99,0.99); 
        RDSP;
    }
}

  /* Motion on a triangulation */
static void
motion_mesh(int position_y, int position_x){
    // mouse motion for rotation
    rotate_x(rx=(float)(position_x-currentx)/Height);
    anglex+=rx;
    rotate_y(ry=(float)(position_y-currenty)/Height);
    angley+=ry;
    currenty=position_y;
    currentx=position_x;
    RDSP;
}

  /* Key input for triangulation */
static void
key_mesh(unsigned char key, int x, int y){

    switch (key){
        case 'd': translate(x0=0.05,0.0,0.0); xx+=x0; RDSP; break;
        case 'a': translate(x0=-0.05,0.0,0.0); xx+=x0; RDSP; break;
        case 'w': translate(0.0,y0=0.05,0.0); yy+=y0; RDSP; break;
        case 's': translate(0.0,y0=-0.05,0.0); yy+=y0; RDSP; break;
        case 'x': translate(0.0,0.0,z0=0.05); zz+=z0; RDSP; break; 
        case 'y': translate(0.0,0.0,z0=-0.05); zz+=z0; RDSP; break; 
        case 0x02B: glScalef(1.01,1.01,1.01); zoom+=0.01; RDSP; break; //'x' key
        case 0x02D: glScalef(0.99,0.99,0.99); zoom-=0.01; RDSP; break; // '-' key
        case 'r': rotate_x(rx=M_PI*0.01); anglex+=rx; RDSP; break;
        case 'f': rotate_x(rx=-M_PI*0.01); anglex+=rx; RDSP; break;
        case 'q': rotate_y(ry=M_PI*0.01); angley+=ry; RDSP; break;
        case 'e': rotate_y(ry=-M_PI*0.01); angley+=ry; RDSP; break;
        case 'h': Printhelp(); break;
        case 'c': glLoadIdentity(); RDSP; break;
        case 'm': meshflag = (meshflag == 'b') ? 'm' : 'b'; RDSP; break;
        case 0x01B: exit(EXIT_SUCCESS); // 'esc'-key
    }
}

/* last but not least the main function */   
int
main(int argc,char **argv){
//   psurface3d sur;
//   int i;             // no idea, what these variables are meant to do
  
  /* Reading mesh */
  if(argc > 1){
    sur0 = read_surface3d(argv[1]); 
    if(argc > 2) meshflag = argv[2][0]; //rudimentary input function to denote, whether to draw mesh or solid
    else meshflag = 'b'; // set solid as standard case
	}
  else{
    printf("No input file!\n");
    return EXIT_SUCCESS;
  }
  if(meshflag != 'm' && meshflag != 'b'){
      Printhelp();
  }else{
      printf("For help press 'h'\n");
  }
    glutInit(&argc, argv);
    glutCreateWindow("Triangulation (for help press 'h')");
    glutPositionWindow(450, 400);
    glutReshapeWindow(800, 800);
	
    glutReshapeFunc(reshape_mesh);
    glutDisplayFunc(display_mesh);

    glutMouseFunc(mouse_mesh);

    glutMotionFunc(motion_mesh);

    glutKeyboardFunc(key_mesh);
    glutMainLoop();
  return EXIT_SUCCESS;
}