/* This program implements an elastic circle composed out of corners which are connected
 * with dampened springs.  An internal pressure is simulated with a simple gamma-gaslaw.
 * 
 * created on 22/01/2003
 * by Jimmy Loodts
 */

#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <list>

#include "Main.h"

// configuration of the corners and the floor
// ------------------------------------------
int		nrCorners		= 6;	// number of corners
double	radiusCircle	= 0.7;	// radius of the circle on which the corners lie
double	floorHeight		= -0.7;	// the height of the floor (0,0) is the center of the circle

// settings for the dynamics
// -------------------------
double	gamma			= 3;	// the power for the gas law
double	pressure		= 2;	// pressure inside the circle
double	cornerMass		= 0.5;	// mass of one corner
double	k				= 200;	// linear spring constant of the contact
double	c				= 70;	// damping constant for the spring
double	g				= -9.81;// gravity acceleration in y-direction

// state variable for the simulation
// ---------------------------------
bool	freeze			= false;	// in order to freeze a certain configuration
bool	stop			= false;	// to stop the simulation
double	time			= 0;		// the simulation time
double	dt				= 0.0008;	// timestep

// the space that is displayed in the window
// -----------------------------------------
double	xmin			= -1.5;
double	xmax			=  1.5;
double	ymin			= -1.5;
double	ymax			=  1.5;

// global declaration of the container
// -----------------------------------
CContainer container(nrCorners, radiusCircle);

// initialize the openGL engine and the container
// ----------------------------------------------
void init(void) 
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glShadeModel(GL_FLAT);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(POINT_SIZE);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth (LINE_WIDTH);
	glEnable(GL_POLYGON_SMOOTH);
	glClear (GL_COLOR_BUFFER_BIT);

	container.initialize(cornerMass, pressure);
}

// draw the scene and swap buffers
// -------------------------------
void display(void)
{
	container.draw();

	glutSwapBuffers();
	glClear (GL_COLOR_BUFFER_BIT);
}

// setting the window frame
// ------------------------
void reshape (int width, int height)
{
   glViewport(0, 0, width, height); 
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluOrtho2D(xmin, xmax, ymin, ymax);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

// calculate the forces
// --------------------
void calculateForces()
{
	container.calculateForces(k, c);
}

// move the container and increase time
// ------------------------------------
void move()
{
	container.move(dt);
	time += dt;
}

// proceed the simulation
// ----------------------
void doCalculations()
{
	if(stop) {}
	else {
		calculateForces();
		if (freeze) container.freeze();
		freeze = false;
		move();
		display();
	}
}

// get user imput fromt the keyboard
// ---------------------------------
void keyboard (unsigned char key, int x, int y)
{
   switch (key) {
		case '+':
			stop = false;
			break;
		case '-':
			stop = true;
			break;
		case 't':
			dt *= 0.5;
			break;
		case 'T':
			dt *= 2;
			break;
		case 'k':
			k *= 0.5;
			break;
		case 'K':
			k *= 2;
			break;
		case 'd':
			c *= 0.9;
			break;
		case 'D':
			c *= 1.1;
			break;
		case 'f':
			freeze = true;
			break;
		case 'p':
			pressure *= 0.5;
			break;
		case 'P':
			pressure *= 2;
			break;
		case 'w':
		//	bool debug = true;
			break;
		case 27: // the escape button
			exit(0);
			break;
		default:
			break;
   }
}

// the main function
// -----------------
int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
   glutInitWindowSize (500, 500); 
   glutInitWindowPosition (100, 100);
   glutCreateWindow (argv[0]);
   init ();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutIdleFunc(doCalculations);
   glutMainLoop();
   return 0;
}