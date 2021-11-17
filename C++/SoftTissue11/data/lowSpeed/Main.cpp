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
double	pressure		= 4;	// pressure inside the circle
double	mass			= 2.85e-9;//0.1;	// mass of one corner
double	k				= 75000.0;//1000;//2000;	// linear spring constant of the contact
double	kCorner			= 0.0;//1000;//2000;	
double	kPressure		= 5e6;//1e4;
double	c				= 2*sqrt(k*mass);//0.004959;//140;	// damping constant for the spring
double	g				= 0.0;//-9.81;// gravity acceleration in y-direction
//double	pullForceX		= 0;	// force (N) with which is pulled in the X-direction
//double	pullForceY		= 0;	// force (N) with which is pulled in the Y-direction
double	displacement		= 0.0;
double	maxStrain		= 5;
double	displacementSpeed	= 5.6e-4;	// m/s increment speed in the position of the northface
double	strain			= 0.0;  // initially there is no strain of the tissue
double	stress				= 0.0;	// the stress measured at the north face
double	loadMax				= 0.0001;	// the maximal load on a cornerPoint at that time
double	volume;
double	initialVolume;

// state variable for the simulation
// ---------------------------------
bool	freeze			= true;	// in order to freeze a certain configuration
bool	stop			= false;	// to stop the simulation
double	time			= 0;		// the simulation time
double	tend			= 0.01*0.01*maxStrain/(double)displacementSpeed; // end time (s) 
						//where 0.01 m is length of the strip
double	dt				= 4e-9;		// timestep
bool	showForces		= false;
int		showTime		= 10000;		// do showTime calculations then display the result
int		nrCalculations	= 0;
int		stretch			= 1;		// start pulling immediately
double	scaleFactor		= 0.1;
int		activeCell		= 0;

// the space that is displayed in the window
// -----------------------------------------
double zoom = 1;
double	xmin			= -zoom*0.0005;
double	xmax			=  zoom*0.01;
double	ymin			= -zoom*0.00525;
double	ymax			=  zoom*0.00525;

// global declaration of the container
// -----------------------------------
CContainer container(mass, pressure);

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
	container.continuePreviousCalculation();

	initialVolume = container.getVolume();
	TimeDelta();
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

	container.calculateForces(time, k, c, displacement, strain, stress);

}

// move the container and increase time
// ------------------------------------
void move()
{
	container.move(dt);

}

// proceed the simulation
// ----------------------
bool stretching = true;
int nrLoops = 0;
double a = (tend/dt)/1000;
int saveInterval = (int)(ceil(a)); // 1000 samples per experiment
void doCalculations()
{

	if(stop) {}
	else {
		if(stretching) {
			if(strain <= maxStrain)	{displacement = displacementSpeed*dt;}
			else {stretching = false;}
		} else {
			if(strain > 0.0) {displacement = -displacementSpeed*dt;}
			else {stretching = true;}
		}
/*	for(kPressure = 1e3; kPressure <1e4; kPressure +=1e2)
	{		
		calculateForces();
		move();
	}// end for(mass = 0.1; mass >0.001; mass -=0.01)
*/
		calculateForces();
	
		// save
		if(nrLoops++ % saveInterval == 0)
		{	
//			cout<<time<<" sec"<<endl;
//			cout.flush();
			const char* nameOut = "data/forceDeformation.txt";
			ofstream outFile;
			outFile.open(nameOut, ios::app);
			if (! outFile){}
			volume = container.getVolume();
			outFile	<<time<<" "
					<<strain<<" "
					<<stress<<" "
					<<100.0*(volume-initialVolume)/initialVolume<<endl;
			outFile.close();
		}

		if (freeze) container.freeze();
		freeze = false;
		if(nrCalculations == showTime){ 
			display();
			nrCalculations = 0;
			cout<<"One step takes "<<TimeDelta()/(double)showTime<<" seconds"<<endl;
			cout<<"Time "<<time<<endl;
		} else {nrCalculations++;}
		move();
		time += dt;
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
		/*case 'f':
			freeze = true;
			cout<<"frozen"<<endl;
			break;
		*/case 'p':
			kPressure *= 0.5;
			cout<<"kPressure: "<<kPressure<<endl;
			break;
		case 'P':
			kPressure *= 2.0;
			cout<<"kPressure: "<<kPressure<<endl;
			break;
		case 'S':
			if(showForces){showForces = false;} else {showForces = true;}
			break;
		case 'f':
			scaleFactor *= 0.5;
			break;
		case 'F':
			scaleFactor *= 2.0;
			break;
		case 'a':
			activeCell--;
			display();
			break;
		case 'A':
			activeCell++;
			display();
			break;
		case 's':
			container.save();
			cout<<"a new inputFile was saved"<<endl;
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
   glutInitWindowSize (900, 500); 
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