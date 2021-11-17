// The debugger can't handle symbols more than 255 characters long.
// STL often creates symbols longer than that.
// When symbols are longer than 255 characters, the warning is disabled.
 #pragma warning(disable:4786)
 #pragma warning(disable:4715) // not all controll paths return a value

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
double	floorHeight		= 0.07;	// the height of the floor (0,0) is the center of the circle

// settings for the dynamics
// -------------------------
double	gamma			= 1.01;	// the power for the gas law
double	pressure		= 1.1;	// pressure inside the circle
double	cornerMass		= 6.147e-9;//0.1;	// mass of one corner
double	k				= 16.66;//1000;//2000;	// linear spring constant of the contact
double	kCorner			= 16.66;//1000;//2000;	
double	c				= 6.4e-4;//0.004959;//140;	// damping constant for the spring
double	g				= 0.0;//-9.81;// gravity acceleration in y-direction
//double	pullForceX		= 0;	// force (N) with which is pulled in the X-direction
//double	pullForceY		= 0;	// force (N) with which is pulled in the Y-direction
double	displacement		= 0.0;
double	displacementSpeed	= 3;	// m/s increment speed in the position of the northface
double	elongation			= 0.0;  // initially there is no elongation of the tissue
double	force				= 0.0;	// the force measured at the north face
double	loadMax				= 0.0001;	// the maximal load on a cornerPoint at that time

// state variable for the simulation
// ---------------------------------
bool	freeze			= false;	// in order to freeze a certain configuration
bool	stop			= true;	// to stop the simulation
double	time			= 0;		// the simulation time
double	tend			= 0.1;		// end time (s)
double	dt				= 1e-6;		// timestep
bool	showForces		= false;
int		stretch			= 1;		// start pulling immediately

// the space that is displayed in the window
// -----------------------------------------
double	xmin			= -0.0025;
double	xmax			=  0.0075;
double	ymin			= -0.002;
double	ymax			=  0.008;

// global declaration of the container
// -----------------------------------
CContainer container(cornerMass, pressure);

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

	container.freeze();
}

// draw the scene and swap buffers
// -------------------------------
void display(void)
{
	container.draw();
	glutSwapBuffers();
	glClear (GL_COLOR_BUFFER_BIT);
	loadMax = 0.0;
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
	container.calculateForces(time, k, c, displacement, elongation, force);

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
bool stretching = true;
int nrLoops = 0;
int saveInterval = (int)((tend/dt)/1000); // 1000 samples per experiment
void doCalculations()
{
	if(stop) {}
	else {
		if(stretching) {
			if(elongation <= 0.003)	{displacement = displacementSpeed*dt;}
			else {stretching = false;}
		} else {
			if(elongation > 0.0006) {displacement = -displacementSpeed*dt;}
			else {stretching = true;}
		}
		
		calculateForces();
		
	
		// save
		if(nrLoops++ % saveInterval == 0)
		{	
			cout<<time<<" sec"<<endl;
			cout.flush();
			const char* nameOut = "data/forceDeformation.txt";
			ofstream outFile;
			outFile.open(nameOut, ios::app);
			if (! outFile){}
			outFile	<<time<<" "
					<<elongation<<" "
					<<force<<endl;
			outFile.close();
		}

		if (freeze) container.freeze();
		freeze = false;
		
		display();move();
	}
}

// get user imput fromt the keyboard
// ---------------------------------
void keyboard (unsigned char key, int x, int y)
{
   switch (key) {
		case '+':
			doCalculations();
			stop = false;
			break;
		case '-':
			stop = true;
			break;
		/*case '4':
			pullForceX -= 5;
			cout<<"pullForceX: "<<pullForceX<<endl;
			break;
		case '6':
			pullForceX += 5;
			cout<<"pullForceX: "<<pullForceX<<endl;
			break;
		*/case '2':
			stretch = -1;
			//pullForceY -= 5;
			//cout<<"pullForceY: "<<pullForceY<<endl;
			break;
		case '8':
			stretch = 1;
			break;
		case '5':
			stretch = 0;
			break;
		case 't':
			dt *= 0.5;
			cout<<"timestep: "<<dt<<endl;
			break;
		case 'T':
			dt *= 2;
			cout<<"timestep: "<<dt<<endl;
			break;
		case 'k':
			k *= 0.5;
			cout<<"spring constant: "<<k<<endl;
			break;
		case 'K':
			k *= 2;
			cout<<"spring constant: "<<k<<endl;
			break;
		case 'd':
			c *= 0.9;
			cout<<"damping constant: "<<c<<endl;
			break;
		case 'D':
			c *= 1.1;
			cout<<"damping constant: "<<c<<endl;
			break;
		case 'f':
			freeze = true;
			cout<<"frozen"<<endl;
			break;
		case 'p':
			gamma *= 0.5;
			cout<<"gamma: "<<gamma<<endl;
			break;
		case 'P':
			gamma *= 2;
			cout<<"gamma: "<<gamma<<endl;
			break;
		case 's':
			showForces = false;
			break;
		case 'S':
			showForces = true;
			break;
/*		case 's':
			container.save();
			break;
*/		case 27: // the escape button
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