// The debugger can't handle symbols more than 255 characters long.
// STL often creates symbols longer than that.
// When symbols are longer than 255 characters, the warning is disabled.
 #pragma warning(disable:4786)
 #pragma warning(disable:985)
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


// settings
// --------
double	g;				// gravity acceleration in y-direction
double	dt;				// timestep
int		forceMode;		//  = 1: pressure + spring
						// != 1: spring
double	time;			// the simulation time in s

bool	proceed;	// to start another simulation
double	scaleFactorForForces; // forces are scaled before drawn with showForces = 1;
bool	showForces;	
int		showTime;		// display and save the results after every showTime steps
int		relaxationTime;	// nr of timesteps for the tissue to relax inbetween to moves
int		timeSteps;	// counter to see when it is showtime
int		nrMoves;	// counter to see the number of times the tissue has been stretched

// the space that is displayed in the window
// -----------------------------------------
double	xmin;
double	xmax;
double	ymin;		
double	ymax;


void readParameterFile(const char* nameIn)
{
	proceed =  true;
	timeSteps = 0;
	nrMoves = 0;
	InputFile iF(nameIn);
	// initialize all the cells with the information strored in the inputfile
	iF.setScope("-parameters-");
	g			= iF.getDouble("g");
	dt			= iF.getDouble("dt");
	forceMode	= iF.getInt("forcemode");
	showTime	= iF.getInt("showtime");
	relaxationTime	= iF.getInt("relaxationtime");
	
	time		= 0;

	scaleFactorForForces	= iF.getDouble("scalefactorforforces");
	showForces	= iF.getBool("showforces");

	// these are set in Container::Container() and BoundaryConditions::Initialize()

/*	xmin		= -0.001; 
	xmax		=  0.018; 
	ymin		=  -0.001; 
	ymax		=  0.018; */
}


// initialize the openGL engine and the parameters
// -----------------------------------------------
//auto//
/*void init() 
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glShadeModel(GL_FLAT);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(POINT_SIZE);
	glEnable(GL_LINE_SMOOTH);
	glLineWidth (LINE_WIDTH);
	glEnable(GL_POLYGON_SMOOTH);
	glClear (GL_COLOR_BUFFER_BIT);

	readParameterFile();

	TimeDelta();
}*/

// draw the scene and swap buffers
// -------------------------------
//auto//
/*void display()
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
}*/

// proceed the simulation
// ----------------------
void doCalculations(Container &container)
{
	container.calculateForces();
	
	if(nrMoves%showTime==0)
	{
//		cout<<"Calculation time per timestep: "<<TimeDelta()/(double)(showTime)<<endl;
		container.saveData(time);
//auto//		display();
//auto//		TimeDelta();
		nrMoves++;
	}
	
	if(timeSteps%relaxationTime==0) {
		if(container.applyBC(1)){proceed = true;}
		else{proceed = false;}
		nrMoves++;
	} else {
		if(container.applyBC(0)){proceed = true;}
		else{proceed = false;};
	}
	container.move(dt);
	
	timeSteps++;

	time += dt;
	
}

// get user imput fromt the keyboard
// ---------------------------------
//auto//
/*void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
		case 'F':
			scaleFactorForForces *= 2;
			break;
		case 'f':
			scaleFactorForForces *= 0.5;
			break;
		case 'R':
			relaxationTime += 1;
			cout<<"relaxationTime = "<<relaxationTime<<endl;
			break;
		case 'r':
			relaxationTime -= 1;
			cout<<"relaxationTime = "<<relaxationTime<<endl;
			break;
		case 's':
			container.save();
			break;
		case 27: // the escape button
			exit(0);
			break;
		default:
			break;
   }
}*/

void getFullName(std::string & fullName, const char* halfName, int number, const char* extension = ".txt")
{
	char temp[33];
	std::string name(halfName);
//	if (number<100) name.append("0");
	if (number<10) name.append("0");
	itoa(number,temp,10);
	name.append(temp);
	name.append(extension);
	fullName = name;
}

void doTestInFolder(int i)
{
	std::string paramStr;
	std::string inputStr;
	std::string outputStr;
	getFullName(paramStr,"data/",i,"/parameterFile.txt");
	getFullName(inputStr,"data/",i,"/long/inputFile.txt");
	getFullName(outputStr,"data/",i,"/long/forceDeformation.txt");
	const char* param  = paramStr.c_str();
	const char* input  = inputStr.c_str();
	const char* output = outputStr.c_str();
	Container container = Container(input, output, param);
   	readParameterFile(param);
	TimeDelta();
    while(proceed){doCalculations(container);}
    cout<<"Time to compute the longitudinal experiment: "<<TimeDelta()<<endl;

	getFullName(inputStr,"data/",i,"/trans/inputFile.txt");
	getFullName(outputStr,"data/",i,"/trans/forceDeformation.txt");
	input  = inputStr.c_str();
	output = outputStr.c_str();
	container = Container(input, output, param);
   	readParameterFile(param);
	TimeDelta();
    while(proceed){doCalculations(container);}
    cout<<"Time to compute the transversal experiment: "<<TimeDelta()<<endl;
}

// the main function
// -----------------
int main(int argc, char** argv)
{
//auto//
	/*glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
   glutInitWindowSize (700, 700); 
   glutInitWindowPosition (300,10);
   glutCreateWindow (argv[0]);
   init ();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutIdleFunc(doCalculations);
   glutMainLoop();
   */
	for(int i=1; i<=6; ++i)
	{
		doTestInFolder(i);
	}
   return 0;
}