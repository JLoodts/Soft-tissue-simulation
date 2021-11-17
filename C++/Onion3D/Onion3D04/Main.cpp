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

double	scaleFactorForForces; // forces are scaled before drawn with showForces = 1;
bool	showForces;	
int		showTime;		// display and save the results after every showTime steps
int		relaxationTime;	// nr of timesteps for the tissue to relax inbetween to moves
int		timeSteps = 0;	// counter to see when it is showtime
int		nrMoves = 0;	// counter to see the number of times the tissue has been stretched

int activeTriangle = 0;

// the space that is displayed in the window
// -----------------------------------------
CTriple center(0.0015,0.001,0); 
CTriple center_begin = center;
CTriple eye(0,0.0001,0.01);
CTriple eye_begin = eye;
CTriple up(0,1,0);
GLfloat light_position[] = {-0.1, 1.2, 0, 1.0};
GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
GLfloat light_diffuse[] = {1, 1, 1, 1.0};
GLfloat light_specular[] = {1, 1, 1, 1.0};
GLfloat spot_direction[] = {1, -1, -0.5};
CPair	spin = ZEROPAIR;
double step = 0.001;
int old_x, old_y; // position of the mouse

double	xmin;
double	xmax;
double	ymin;		
double	ymax;


void readParameterFile()
{
	const char* nameIn = "data/parameterFile.txt";
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

}

// global declaration of the container
// -----------------------------------
Container container = Container();

// initialize the openGL engine and the parameters
// -----------------------------------------------
void init() 
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
}

// draw the scene and swap buffers
// -------------------------------
void display()
{
/*	container.draw();
	glutSwapBuffers();
	glClear (GL_COLOR_BUFFER_BIT);
	*/
	// eye = Bezier(eye0,eye1,eye2,eye3,(float)time/100);
	

		
//		glColor4f(0,1,0,0.2);
//		drawPlane(CTriple(10,0,0), CTriple(0,0,10), CTriple(-5,0,-5));

	glPopMatrix();

	glutSwapBuffers();	
	
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	glPushMatrix();
		gluLookAt(eye.x,eye.y,eye.z,center.x,center.y,center.z,up.x,up.y,up.z);
		glRotatef(spin.y, 1, 0, 0);
		glRotatef(spin.x, 0, 1, 0);
		container.draw();
}

// setting the window frame
// ------------------------
void reshape (int width, int height)
{

	  float black[] = { 0, 0, 0, 0 };

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30, (float)width/height, 0.001, 1000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity(); 
    glFogfv(GL_FOG_COLOR, black);
    glFogf(GL_FOG_START, 5);
    glFogf(GL_FOG_END, 6);
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, GL_LINEAR);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(POINT_SIZE);
	//glEnable(GL_LINE_SMOOTH);
	glLineWidth (LINE_WIDTH);
	//glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_DEPTH_TEST);
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient); 
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse); 
	//glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular); 
	glLightfv(GL_LIGHT0, GL_POSITION, light_position); 
	//glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.5);
	//glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.5);
	//glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.2);
	//glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 89.0);
	//glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spot_direction);
	//glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 2.0);
	glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
}



// proceed the simulation
// ----------------------
void doCalculations()
{
	container.calculateForces();
	
	if(nrMoves%showTime==0)
	{
//		cout<<"Calculation time per timestep: "<<TimeDelta()/(double)(showTime)<<endl;
		container.saveData(time);
		display();
		TimeDelta();
		nrMoves++;
	}
	
	if(timeSteps%relaxationTime==0) {
		container.applyBC(1);
		nrMoves++;
	} else {
		container.applyBC(0);
	}

	container.move(dt);
	
	timeSteps++;

	time += dt;
	
}

// get user imput fromt the keyboard
// ---------------------------------
void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
	    case 'a':
			if(activeTriangle>0){activeTriangle--;}
			break;
		case 'A':
			if(activeTriangle<100){activeTriangle++;}
			break;
		case 'F':
			scaleFactorForForces *= 2;
			cout<<"scaleFactorForce = "<<scaleFactorForForces<<endl;
			break;
		case 'f':
			scaleFactorForForces *= 0.5;
			cout<<"scaleFactorForce = "<<scaleFactorForForces<<endl;
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
		case '8': 
			eye.y+=step; 
			center.y+=step;;         
			break;
		case '2': 
			eye.y-=step; 
			center.y-=step;         
			break;
		case '6': 
			eye.x+=step;
			center.x+=step;               
			break;
		case '4': 
			eye.x-=step; 
			center.x-=step;         
			break;
		case '1': 
			eye.z+=step; 
			center.z+=step;
			break;
		case '7': 
			eye.z-=step;
			center.z-=step;
			break;
		case '*': 
			center.x+=step;                 
			break;
		case '/': 
			center.x-=step;    
			break;
		case '-': 
			center.y+=step;               
			break;
		case '+': 
			center.y-=step;           
			break;
		case '3': 
			center.z+=step;               
			break;
		case '9': 
			center.z-=step;           
			break;
		case '5': 
			eye = eye_begin; 
			center = center_begin;
			up = CTriple(0,1,0);          
			break;
		case 'g':
			step *= 2;
			break;
		case 'd':
			step /=2;
			break;
		case 27: // the escape button
			exit(0);
			break;
		default:
			break;
   }
}

void mouse(int button, int state, int x, int y)
{
    old_x = x;
    old_y = y;

    glutPostRedisplay();
}

void motion(int x, int y)
{
    spin.x = x - old_x;
    spin.y = y - old_y;

    glutPostRedisplay();
}

// the main function
// -----------------
int main(int argc, char** argv)
{
   glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_DOUBLE);
   glutInitWindowSize (700, 700); 
   glutInitWindowPosition (300,10);
   glutCreateWindow (argv[0]);
   init ();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutMotionFunc(motion);
   glutMouseFunc(mouse);
   glutKeyboardFunc(keyboard);
   glutIdleFunc(doCalculations);
   glutMainLoop();
   return 0;
}