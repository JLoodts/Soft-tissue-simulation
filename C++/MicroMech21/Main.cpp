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
#include "General.h"
#include "Draw.h"
#include "Ball.h"
#include "Container.h"

// configuration of the balls and the walls
// ----------------------------------------
int		maxNrBalls		= 10000;		// number of balls
double	minRadius		= 0.0003;	// minimum radius of the balls
double	maxRadius		= 0.0004;// radius of the balls

/*
 *
 *             Y ^
 *               |
 *        ceilY -|          -----------------------------------
 *               |   INLET |->                                 |
 *  inletFloorY _|         |______________                     | OUTLET
 *               |                        |                    |
 * outletFloorY _|                        |____________________|
 *               |
 *                ---------|--------------|--------------------|------>
 *                      inletX         tresholdX           outletX    X
 *
 */
double	ceilY			= 0.01;	// upper point in y direction
double	inletFloorY		= 0.005;	// lower point in y direction at the inlet
double	outletFloorY	= 0.0;		// lower point in y direction at the outlet
//double	amplitude		= -0.004;
//double	frequency		= 30.0;
double	inletX			= 0.0;		// inlet in x direction
double	tresholdX		= 0.02;// treshold position in x direction
double	outletX			= 0.08;// outlet in x direction


// settings for the dynamics
// -------------------------
double	density			= 1.23;		// kg/m³
double	minMass			= 5e-6;
double	maxMass			= 5e-6;
//double	ballMass		= 1e-4;//3.5*ballRadius*ballRadius*density;		// mass of one ball
double	k				= 2000;	//500 for water linear spring constant of the contact
double	c				= 0.30;//0.07 for water 2.0*sqrt(maxMass*k);		// damping constant for the spring
double	vMean			= 0.1;	// mean velocity in m/s
double	angle			= 180;		// angle in which gravity works in degrees
CPair	gravity			= ZEROPAIR;//CPair(0.0,-GRAVITY);
double	swarmImpact		= 0.0007;		// determines how big the influence of the swarm is
double friction			= 0.9;			// % of vel that remains after collision

// state variable for the simulation
// ---------------------------------
bool	stop			= true;	// autorun on or off
int		nrShows			= 0;		// to organise when we have to draw the scene
const int showTime		= 40;		// draw every 10 calculations
double	time			= 0;		// the simulation time
double	dt				= 1e-3;	// timestep
double	dt_impact		= (PI/sqrt(k/minMass-c/minMass))/10.0;	// 10 steps for one collision function of m and k
double	dt_overlap;					// adjust timestep to max speed
int timeSteps = 0;
double elapsedTime = 0;
// the space that is displayed in the window
// -----------------------------------------
double	xmin			= inletX-0.005;
double	xmax			= outletX+0.005;
double	ymin			= outletFloorY-0.005;//0.5*(outletX-inletX);
double	ymax			= ceilY+0.005;//0.5*(outletX-inletX);
const char* fileName = "output";

// global declaration of the container
// -----------------------------------
CContainer container(maxNrBalls, minRadius, maxRadius, minMass, maxMass);

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

	container.initialize();

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
double velMax = 0;
void move()
{
	
	velMax = container.move(dt);
	double	impactDuration	= PI/sqrt(k/minMass-c/minMass);
	dt_impact	= impactDuration/20.0;
	dt_overlap	= 0.07*minRadius/(1.42*velMax);	// *1.42 omdat slechts 1 richting werd beschouwd in velMax
//out<<"dt_overlap: "<<dt_overlap<<" dt_impact: "<<dt_impact<<endl;
//	cout<<"max vel: "<<velMax<<"  dt: "<<dt<<endl;
	if(dt_overlap<dt_impact){dt = dt_overlap;}else{dt = dt_impact;}
	time += dt;
//	outletFloorY = amplitude*cos(frequency*time*2*PI);
//	glColor3f(1,0,0); DrawLine(inletX,outletFloorY+0.001, outletX, outletFloorY+0.001);
}

// proceed the simulation
// ----------------------
double timeCounts = 0;
double sum = 0;
void doCalculations()
{
	if(stop) {}
	else {
		TimeDelta();
		calculateForces();
		move();
		elapsedTime += TimeDelta();
		timeSteps++;
		/*if(timeSteps>=1)
		{
			timeCounts++;
			double compTime = elapsedTime/(double)(timeSteps);
			sum += compTime;
			cout<<"nr balls: "<<container.getNrBalls()<<endl
				<<"computing time: "<<compTime<<endl
				<<"mean comp time: "<<sum/timeCounts<<endl;
			timeSteps = 0;
			elapsedTime = 0;
		}*/
		if(nrShows == showTime){ 
			display();
			nrShows = 0;
		} else {nrShows++;}
	}
}


// get user imput fromt the keyboard
// ---------------------------------
void keyboard (unsigned char key, int x, int y)
{
   switch (key) {
		case '+':
			doCalculations();
			break;
		case 't':
			dt *= 0.5;
			cout<<"Timestep :"<<dt<<endl;
			break;
		case 'T':
			dt *= 2;
			cout<<"Timestep :"<<dt<<endl;
			break;
		case 'k':
			k *= 0.5;
			cout<<"Spring constant :"<<k<<endl;
			break;
		case 'K':
			k *= 2;
			cout<<"Spring constant :"<<k<<endl;
			break;
		case 'd':
			c *= 0.9;
			cout<<"Damping constant :"<<c<<endl;
			break;
		case 'D':
			c *= 1.1;
			cout<<"Damping constant :"<<c<<endl;
			break;
/*		case 'A':
			amplitude *= 1.1;
			cout<<"Amplitude :"<<amplitude<<endl;
			break;
		case 'a':
			amplitude *= 0.9;
			cout<<"Amplitude :"<<amplitude<<endl;
			break;
		case 'F':
			frequency *= 1.1;
			cout<<"Frequency :"<<frequency<<endl;
			break;
		case 'f':
			frequency *= 0.9;
			cout<<"Frequency :"<<frequency<<endl;
			break;
*/		case 'c':
			container.clearVelField();
			cout<<"The velocity field is cleared"<<endl;
			break;
		case 'w':
			swarmImpact *= 0.9;
			cout<<"swarm impact: "<<swarmImpact<<endl;
			break;
		case 'W':
			swarmImpact *= 1.1;
			cout<<"swarm impact: "<<swarmImpact<<endl;
			break;
		case 's':
			container.saveVelField(fileName);
			cout<<"The velocity field has been saved"<<endl;
			break;
		case 'p':
			stop = true;
			cout<<"Program execution is paused press 'r' to resume"<<endl;
			break;
		case 'r':
			stop = false;
			cout<<"Program execution resumed, press 's' to pause"<<endl;
			break;
		case 'i':
			cout<<"Info about the program settings"<<endl
				<<"--------------------------------"<<endl
				<<"angle: "<<angle<<" gX : "<<gravity.x<<"   gY: "<<gravity.y<<endl
				<<"nr balls: "<<container.getNrBalls()<<endl
				<<"nr cells: "<<container.getNrCells()<<endl
				<<"Timestep :"<<dt<<endl
				<<"dt_overlap: "<<dt_overlap<<" dt_impact: "<<dt_impact<<endl
				<<"Spring constant :"<<k<<endl
				<<"Damping constant :"<<c<<endl
				<<"computing time: "<<elapsedTime/(double)(timeSteps)<<endl
				<<"press 's' to pause program execution"<<endl
				<<"press 'r' to resume"<<endl
				<<"--------------------------------"<<endl;
			timeSteps = 0; elapsedTime = 0; timeCounts = 0; sum = 0;
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
   glutInitWindowSize (900,200); 
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