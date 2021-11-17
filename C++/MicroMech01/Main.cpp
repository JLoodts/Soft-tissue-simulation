#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535898

#include "Draw.h"
#include "Contact.h"

GLuint startList;
float **wall;
float *stirrer;
float **disc;
const int w = 75;	// number of corners
const int d = 2;	// number of discs
const float rw = 1;	// radius of wallcircle
const float rd = 0.3;	// radius of a disc
const float md = 1;	// mass of one disc
const float mw = 1;	// mass of one wall
const float stirrSpeed = 100; // speed of the stirrer
float stirrRadius = 0.5;

float time = 0;
float k = 20000;	// spring constant of the contact
float damp = 0.99;		// damping coefficient of the wall-spring
float dt = 0.005;  // timestep
const float g = 0; // gravity acceleration

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

	// create an array of n x 4 elements to store wall positions
	wall = new float*[w];
	for(int i=0; i<w; i++)
	{
		wall[i] = new float[8];
		for(int j=0; j<8; j++) wall[i][j] = 0;
	}

	// create the stirrer
	stirrer = new float[4];
	for(int j=0; j<4; j++) stirrer[j] = 0;

	// create an array of d x 2 elements to store disc positions
	disc = new float*[d];
	for(i=0; i<d; i++)
	{
		disc[i] = new float[6];
		for(int j=0; j<6; j++) disc[i][j] = 0;
	}

	// calculate beginpositions of an n-polygon of walls 
	// around the origin
	float x1, x2, y1, y2;
	x1 = rw; y1 = 0;
	for(i=0; i<w; i++)
	{
		float angle = 2*PI*(i+1)/(float)(w);
		x2 = rw*cos(angle); y2 = rw*sin(angle);
		wall[i][0] = x1; wall[i][1] = y1;
		wall[i][2] = x2; wall[i][3] = y2;
		x1 = x2; y1 = y2;
	}
	
	// calculate beginposition of the stirrer
	stirrer[0] = stirrRadius*sin(stirrSpeed*0.5*time/PI);
	stirrer[1] = stirrRadius*cos(stirrSpeed*0.5*time/PI);
	stirrer[2] = -stirrRadius*sin(stirrSpeed*0.5*time/PI);
	stirrer[3] = -stirrRadius*cos(stirrSpeed*0.5*time/PI);

	// calculate beginpositions of the discs
	disc[0][0] = 0.5; disc[0][1] = 0;
	disc[1][0] = -0.5; disc[1][1] = 0;
	// beginvelocities for the discs
	disc[0][2] = -1; disc[0][3] = -0.2;
	disc[1][2] = 1; disc[1][3] = 0.2;
}

void display(void)
{
	

	// draw the discs in blue
	glColor3f(0,0,1);
	drawDiscs(disc, d, rd);

	// draw all the boundary lines in blue
	glColor3f(0, 0, 1);
	drawWalls(wall, w);

	// draw all the springs between the walls in purple
	glColor3f(1, 0, 1);
	drawSprings(wall, w);

	// draw the stirrer in black
	glColor3f(0, 0, 0);
	glBegin(GL_LINES); 
		glVertex2f(stirrer[0], stirrer[1]); 
		glVertex2f(stirrer[2], stirrer[3]); 
	glEnd();

	glColor3f(1,0,0);
	glutSwapBuffers();
	glClear (GL_COLOR_BUFFER_BIT);
}

void reshape (int width, int height)
{
   glViewport(0, 0, width, height); 
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 20.0);
   //gluPerspective(30, (GLdouble)width/(GLdouble)height, 0.1, 10);
   //glOrtho(-1.5, 1.5, -1.5, 1.5, 1, 10);
   gluOrtho2D(-1.5, 1.5, -1.5, 1.5);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   //gluLookAt(0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);   
}


void Contactdetection()
{
	/*
	 *	All contacts between walls and disc are checked
	 *	if a contact occors the resulting force is taken into 
	 *	account as accelerations of the objects in contact
	 */

	float xclosest, yclosest, xnew, ynew, xw, yw, xwunit, ywunit,
		l, lw, distance, overlap, dirx, diry, forcex, forcey;
	for(int i=0; i<d; ++i)
	{
		for(int j=0; j<w; ++j)
		{
			// find closest point of the wall to the disc
			xnew = disc[i][0] - wall[j][0];
			ynew = disc[i][1] - wall[j][1];
			xw = wall[j][2] - wall[j][0];
			yw = wall[j][3] - wall[j][1];
			lw = sqrt(xw*xw + yw*yw);
			xwunit = xw/lw;
			ywunit = yw/lw;
			l = (xnew*xwunit + ynew*ywunit);
			if(l<=0)
			{
				xclosest = wall[j][0]; yclosest = wall[j][1];
			}
			else
			{
				if(l >= lw)
				{
					xclosest = wall[j][2]; yclosest = wall[j][3];
				}
				else
				{
					xclosest = l*xwunit + wall[j][0];
					yclosest = l*ywunit + wall[j][1];
				}
			}
			

			// calculate the force and modify the accelerations
			dirx = disc[i][0]-xclosest; diry = disc[i][1]-yclosest;
			distance = sqrt(dirx*dirx + diry*diry);
			overlap = distance - rd;
			if(overlap<0)
			{
				drawPoint(xclosest,yclosest);

				forcex = -overlap*k*(dirx)/distance;
				forcey = -overlap*k*(diry)/distance;
				disc[i][4] += forcex;// for not-unity mass/md;
				disc[i][5] += forcey;// for not-unity mass/md;
				wall[j][6] += -forcex;// for not-unity mass/mw;
				wall[j][7] += -forcey;// for not-unity mass/mw;
			}
		}
	}

	// model springs connecting the walls
	for(i=0; i<w-1; ++i)
	{
		// calculate the force and modify the accelerations
		dirx = wall[i+1][0]-wall[i][2]; diry = wall[i+1][1]-wall[i][3];
		forcex = k*dirx;
		forcey = k*diry;
		wall[i][6] += forcex;// for not-unity mass/mw;
		wall[i][7] += forcey;// for not-unity mass/mw;
		wall[i+1][6] += -forcex;// for not-unity mass/mw;
		wall[i+1][7] += -forcey;// for not-unity mass/mw;
	}
	// to connect the last wall back to the first
	dirx = wall[0][0]-wall[w-1][2]; diry = wall[0][1]-wall[w-1][3];
	forcex = k*dirx;
	forcey = k*diry;
	wall[w-1][6] += forcex;// for not-unity mass/mw;
	wall[w-1][7] += forcey;// for not-unity mass/mw;
	wall[0][6] += -forcex;// for not-unity mass/mw;
	wall[0][7] += -forcey;// for not-unity mass/mw;

	// contact forces between discs
	for(i=0; i<d-1; ++i)
	{
		for(int j=i+1; j<d; ++j)
		{
			// calculate the force and modify the accelerations
			dirx = disc[j][0]-disc[i][0]; diry = disc[j][1]-disc[i][1];
			distance = sqrt(dirx*dirx + diry*diry);
			overlap = distance - 2*rd;
			if(overlap<0)
			{
				drawPoint(disc[i][0]+0.5*dirx,disc[i][1]+0.5*diry);

				forcex = -overlap*k*(dirx)/distance;
				forcey = -overlap*k*(diry)/distance;
				disc[i][4] += -forcex;// for not-unity mass/mw;
				disc[i][5] += -forcey;// for not-unity mass/mw;
				disc[j][4] += forcex;// for not-unity mass/mw;
				disc[j][5] += forcey;// for not-unity mass/mw;
			}
		}
	}

	// contact forces between discs and stirrer
	for(i=0; i<d; ++i)
	{
		// find closest point of the stirrer to the disc
		xnew = disc[i][0] - stirrer[0];
		ynew = disc[i][1] - stirrer[1];
		xw = stirrer[2] - stirrer[0];
		yw = stirrer[3] - stirrer[1];
		lw = sqrt(xw*xw + yw*yw);
		xwunit = xw/lw;
		ywunit = yw/lw;
		l = (xnew*xwunit + ynew*ywunit);
		if(l<=0)
		{
			xclosest = stirrer[0]; yclosest = stirrer[1];
		}
		else
		{
			if(l >= lw)
			{
				xclosest = stirrer[2]; yclosest = stirrer[3];
			}
			else
			{
				xclosest = l*xwunit + stirrer[0];
				yclosest = l*ywunit + stirrer[1];
			}
		}
		
		// calculate the force and modify the accelerations
		dirx = disc[i][0]-xclosest; diry = disc[i][1]-yclosest;
		distance = sqrt(dirx*dirx + diry*diry);
		overlap = distance - rd;
		if(overlap<0)
		{
			drawPoint(xclosest,yclosest);

			forcex = -overlap*k*(dirx)/distance;
			forcey = -overlap*k*(diry)/distance;
			disc[i][4] += forcex;// for not-unity mass/md;
			disc[i][5] += forcey;// for not-unity mass/md;
		}
	}
}

void Move()
{
	for(int i=0; i<d; ++i)
	{
		disc[i][2] += disc[i][4]*dt;
		disc[i][3] += (disc[i][5] - g)*dt;
		disc[i][2] *= damp;
		disc[i][3] *= damp;
		disc[i][0] += disc[i][2]*dt;
		disc[i][1] += disc[i][3]*dt;
		disc[i][4] = 0;
		disc[i][5] = 0;
	}

	// keep wall nr 0 fixed
	wall[0][6] = 0; wall[0][7] = 0;
	wall[0][4] = 0; wall[0][5] = 0;

	for(i=0; i<w; ++i)
	{
		wall[i][4] += wall[i][6]*dt;
		wall[i][5] += wall[i][7]*dt;
		wall[i][4] *= damp;
		wall[i][5] *= damp;
		wall[i][0] += wall[i][4]*dt;
		wall[i][1] += wall[i][5]*dt;
		wall[i][2] += wall[i][4]*dt;
		wall[i][3] += wall[i][5]*dt;
		wall[i][6] = 0;
		wall[i][7] = 0;
	}

	// rotate the stirrer
	stirrer[0] = stirrRadius*sin(stirrSpeed*0.5*time/PI);
	stirrer[1] = stirrRadius*cos(stirrSpeed*0.5*time/PI);
	stirrer[2] = -stirrRadius*sin(stirrSpeed*0.5*time/PI);
	stirrer[3] = -stirrRadius*cos(stirrSpeed*0.5*time/PI);

	time += dt;
}

void DoCalculations()
{
	display();
	Contactdetection();
	//glutPostRedisplay();
	Move();
}

void keyboard (unsigned char key, int x, int y)
{
   switch (key) {
		case '+':
			DoCalculations();
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
			damp *= 0.99;
			break;
		case 'D':
			damp *= 1.01;
			break;
		case 27:
			exit(0);
			break;
		default:
			break;
   }
}

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
   //glutIdleFunc(DoCalculations);
   glutMainLoop();
   return 0;
}
