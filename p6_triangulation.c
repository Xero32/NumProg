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
#define RDSP glutPostRedisplay();
#ifndef M_PI
/* Copied from math.h  */
#define M_PI		3.14159265358979323846
#endif
  
 

/* ------------------------------------------------------------
* Global variable
*-------------------------------------------------------------*/
							
psurface3d sur0;//sur1; 
real x0 = -0.0, y0 = 0.0, z0 = -0.0; //coordinates of mesh
// real x1 = 2.0, y1 = 2.0, z1 = -10.0; //coordinates of solid
real xx = 0.0, yy = 0.0, zz = -5.0;
// real xtr,ytr,ztr; // x-translation, y-translation, z-translation aka zoom
real rx,ry, anglex, angley;
int k = 0;
int currentx, currenty;
int Height, Width;
int rotateflag;
char meshflag;
double zoom = 1.0;
/* Translation */

static void
Printhelp(){
    printf("\nUsage Info:\n");
    printf("After inputfile, denote simulation mode:\n\t'm': mesh\n\t'b': solid\n\n");
    printf("Use Keyboard to move, zoom, and rotate objects.\n\n");
    printf("Use 'w', 'a', 's', 'd' to move object.\nZoom in with '+', zoom out with '-'\n");
    printf("Rotate around the x axis by pressing 'r' and 'f'\n");
    printf("Rotate around the y axis by pressing 'q' and 'e'\n\n");
    printf("You can also rotate the body with your mouse by pressing any mouse button.\nZoom with mouse wheel\n\n");
    printf("Close with 'esc'\n\n");
}


static void
translate(double x, double y, double z){
//     glLoadIdentity();
//     glPushMatrix();
    
    double tr[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, x,y,z,1};
    glMultMatrixd(tr);
//     glutPostRedisplay();
}

/* Rotation around x-axis */
static void
rotate_x(double alpha){
//     glLoadIdentity();
//     glPushMatrix();
    double rotx[16] = {1,0,0,0, 0,cos(alpha),sin(alpha),0, 0,-sin(alpha),cos(alpha),0, 0,0,0,1};
    glMultMatrixd(rotx);
//     glutPostRedisplay();    
}

/* Rotation around y-axis */
static void
rotate_y(double alpha){
//     glLoadIdentity();
//     glPushMatrix();
    double roty[16] = {cos(alpha),0,-sin(alpha),0, 0,1,0,0, sin(alpha),0,cos(alpha),0, 0,0,0,1};
    glMultMatrixd(roty);
//     glutPostRedisplay();
}


  
 /* Drawing complete surface triangulation */
static void
display_mesh(){
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_DEPTH_TEST);
    printf("%f, %f, %f\t\t %f, %f\n",xx,yy,zz,rx,ry);
    printf("%f, %f, %f\t\t %f, %f\n",x0,y0,z0,anglex,angley);
    printf("%f\n",zoom);
    
    glClearColor(0.9,0.9,0.9,1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
    glPushMatrix();
    
    if(rotateflag){
        rotate_y(rx);
        rotate_x(ry);
    }else{
        rotate_x(rx);
        rotate_y(ry);
    }
    translate(x0,y0,z0);
    glPopMatrix();
    printf("counter: %d\n",k++);
    glColor3f(0.3,0.2,0.4);
    if(meshflag == 'm'){
        glBegin(GL_LINES);
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        
        int edges = sur0->edges;
        int (*e)[2] = sur0->e;
        real (*x)[3] = sur0->x;
        
        for(int i = 0; i < edges; i++){
            for(int j = 0; j < 2; j++){
                glVertex3f(x[e[i][j]][0], x[e[i][j]][1], x[e[i][j]][2]);
            }
        }
        
    }else if(meshflag == 'b'){
        int triangles=sur0->triangles;
	int (*t)[3]=sur0->t;
	real (*x)[3]=sur0->x;
	real (*n)[3]=sur0->n;
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
    /*
    if(width > height)
        glViewport((width-height)/2, 0, height, height);
    else
        glViewport(0, (height-width)/2, width, width);
    */
    printf("RESHAPE\n");
    gluPerspective(35.0, 1.0, 1.0, 30.0);
    glTranslatef(xx,yy,zz);
    rotate_y(angley);
    rotate_x(anglex);
    
//     glRotatef(anglex,1.0,0.0,0.0);
//     glRotatef(angley,0.0,1.0,0.0);
//     translate(x0,y0,z0);
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

  static void
mouse_mesh(int button, int state, int position_y, int position_x){
    currentx = position_x;
    currenty = position_y;
    if(button == 3){
        translate(0.0,0.0,z0=0.05);
        zz+=z0;
        RDSP;
    }else if(button == 4){
        translate(0.0,0.0,z0=-0.05);
        zz+=z0;
        RDSP;
    }
}

  /* Motion on a triangulation */
static void
motion_mesh(int position_y, int position_x){
    rotate_x(rx=(float)(position_x-currentx)/Height);
    anglex+=rx;
    rotate_y(ry=(float)(position_y-currenty)/Height);
    angley+=ry;
    currenty=position_y;
    currentx=position_x;
    glutPostRedisplay();
}

 
  /* Key input for triangulation */
static void
key_mesh(unsigned char key, int x, int y){

    switch (key){
        case 'd': translate(x0=0.05,0.0,0.0); xx+=x0; RDSP; break;
        case 'a': translate(x0=-0.05,0.0,0.0); xx+=x0; RDSP; break;
        case 'w': translate(0.0,y0=0.05,0.0); yy+=y0; RDSP; break;
        case 's': translate(0.0,y0=-0.05,0.0); yy+=y0; RDSP; break;
        case 0x02B: translate(0.0,0.0,z0=0.05); zz+=z0; RDSP; break; //'x' key
        case 0x02D: translate(0.0,0.0,z0=-0.05); zz+=z0; RDSP; break; // '-' key
        
        case 'r': rotate_x(rx=M_PI*0.01); anglex+=rx; rotateflag = 0; RDSP; break;
        case 'f': rotate_x(rx=-M_PI*0.01); anglex+=rx; rotateflag = 0; RDSP; break;
        case 'q': rotate_y(ry=M_PI*0.01); angley+=ry; rotateflag = 1; RDSP; break;
        case 'e': rotate_y(ry=-M_PI*0.01); angley+=ry; rotateflag = 1; RDSP; break;
        case 'h': Printhelp(); break;

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
    if(argc > 2) meshflag = argv[2][0];
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

//     glutTimerFunc(50, timer, 0);
	
	glutMainLoop();
  return EXIT_SUCCESS;
}