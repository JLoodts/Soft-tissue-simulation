/* 
 *  Demonstration of picking and rendering luminous objects.  Drag the
 *  middle mouse button to spin the object.  Move the mouse over the
 *  bulbs to light them.
 * 
 *  author: Nate Robins
 *  email: ndr@pobox.com
 *  www: http://www.pobox.com/~ndr 
 */
#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include <GL/glut.h>

#include "Container.h"

double	DISC_SIZE		= 0.1;
Pair	mousePos		(  0,	0  );
GLint	mouse_state		= -1;
GLint	mouse_button	= -1;
GLuint  windowWidth		= 512;				// in pixels
GLuint  windowHeight	= 512;
Pair	dimMin				(  0,	 0  );
Pair	dimMax				( 0.003,	0.003 );
double	tissueThickness	= 120e-6;
Pair	nrGridPoints	( 71,	13 );
int		selectedCornerNr = 0;
const	GLsizei SELECT_BUFFER = 30000;
GLuint	select_buffer[SELECT_BUFFER];	// selection buffer	
int		begin			= 0;
int		end				= 0;
bool	floating		= false;
bool	showPBSprings	= true;
bool	showCBSprings	= true;
bool	showCells		= true;

Container container(nrGridPoints);
void keyboard(unsigned char key, int x, int y);

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
	glSelectBuffer(SELECT_BUFFER, select_buffer); // set the selection buffer size 
	keyboard('h',0,0);
}

void reshape(int width, int height)
{
	glViewport(0, 0, width, height); 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(dimMin.x, dimMax.x, dimMin.y, dimMax.y);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	container.draw();
	if((floating)&&(begin)){
		glColor3f(1,0,1);
		drawLine(Pair(container.getCornerPos(begin)),Pair(mousePos.x,mousePos.y));
	}

	glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'h':
		switch(mode){
		case CBSPRING:
			cout<<"--- Cell Boundary Spring Mode ---"<<endl;
			cout<<"    -------------------------    "<<endl<<endl;
			break;
		case PBSPRING:
			cout<<"--- Plane Boundary Spring Mode ---"<<endl;
			cout<<"    --------------------------    "<<endl<<endl;
			break;
		case CELL:
			cout<<"--- Cell Mode ---"<<endl;
			cout<<"    ---------    "<<endl<<endl;
			break;
		case NORTH_AND_SOUTH:
			cout<<"--- North and South Mode ---"<<endl;
			cout<<"    --------------------    "<<endl<<endl;
			break;
		case EAST:
			cout<<"--- Eastern Boundary Mode ---"<<endl;
			cout<<"    ---------------------    "<<endl<<endl;
			break;
		case WEST:
			cout<<"--- Western Boundary Mode ---"<<endl;
			cout<<"    ---------------------    "<<endl<<endl;
			break;
		case CBTRANS:
			cout<<"--- Cell Boundary Transsection Mode ---"<<endl;
			cout<<"    -------------------------------    "<<endl<<endl;
			break;
		case PBTRANS:
			cout<<"--- Plane Boundary Transsection Mode ---"<<endl;
			cout<<"    --------------------------------    "<<endl<<endl;
			break;
		}
		printf("1            -  Toggle display cell boundary springs\n");
		printf("2            -  Toggle display plane boundary springs\n");
		printf("3            -  Toggle display cells\n");
		printf("s            -  Save to savedInputFile.txt\n");
		printf("e            -  Export as inputFile.txt\n");
		printf("l            -  Load inputFile.txt\n");
		printf("m            -  Mode switch\n");
		printf("u            -  Undo\n");
		printf("backspace    -  Reset\n");
		printf("escape or q  -  Quit\n\n");
		break;
	case 'm':
		switch(mode){
		case CBSPRING:
			mode = PBSPRING;
			cout<<"--- Plane Boundary Spring Mode ---"<<endl;
			cout<<"    --------------------------    "<<endl<<endl;
			break;
		case PBSPRING:
			mode = CELL;
			cout<<"--- Cell Mode ---"<<endl;
			cout<<"    ---------    "<<endl<<endl;
			break;
		case CELL:
			mode = NORTH_AND_SOUTH;
			cout<<"--- North and South Mode ---"<<endl;
			cout<<"    --------------------    "<<endl
				<<"First select the south point (orange) then the north point (blue)."<<endl<<endl;
			break;
		case NORTH_AND_SOUTH:
			mode = EAST;
			cout<<"--- Eastern Boundary Mode ---"<<endl;
			cout<<"    ---------------------    "<<endl<<endl;
			break;
		case EAST:
			mode = WEST;
			cout<<"--- Western Boundary Mode ---"<<endl;
			cout<<"    ---------------------    "<<endl<<endl;
			break;
		case WEST: 
			mode = CBTRANS;
			cout<<"--- Cell Boundary Transsection Mode ---"<<endl;
			cout<<"    -------------------------------    "<<endl<<endl;
			break;
		case CBTRANS:
			mode = PBTRANS;
			cout<<"--- Plane Boundary Transsection Mode ---"<<endl;
			cout<<"    --------------------------------    "<<endl<<endl;
			break;
		case PBTRANS:
			mode = CBSPRING;
			cout<<"--- Cell Boundary Spring Mode ---"<<endl;
			cout<<"    -------------------------    "<<endl<<endl;
			break;
		}
		break;
	case 'u':
		switch(mode){
		case CBSPRING:
			container.removeLastSpring();
			cout<<"The last cell boundary spring is removed"<<endl;
			break;
		case PBSPRING:
			container.removeLastSecSpring();
			cout<<"The last plane boundary spring is removed"<<endl;
			break;
		case CELL:
			container.removeLastCell();
			cout<<"The last cell is removed"<<endl;
			break;
		case NORTH_AND_SOUTH:
			container.resetNorthAndSouth();
			cout<<"North and South are reset"<<endl;
			break;
		case EAST:
			container.removeLastEastCorner();
			cout<<"The last east boundary corner is removed"<<endl;
			break;
		case WEST:
			container.removeLastWestCorner();
			cout<<"The last west boundary corner is removed"<<endl;
			break;
		case CBTRANS:
			container.removeLastStressSpring();
			cout<<"The last cell boundary stress spring is removed"<<endl;
			break;
		case PBTRANS:
			container.removeLastStressSecSpring();
			cout<<"The last plane boundary stress spring is removed"<<endl;
			break;
		}
		break;
	case 'l':
		container.load();
		cout<<"The configuration from 'inputFile.txt' was loaded"<<endl;
		break;
	case 's':
		container.save();
		cout<<"All data were saved to 'inputFile.txt' in the /data directory"<<endl;
		break;
	case 'e':
		container.save4SoftTissue();
		cout<<"All data were saved to 'inputFile4SoftTissue.txt' in the /data directory"<<endl;
		break;
	case '1':
		showCBSprings = !showCBSprings;
		if(showCBSprings){
			cout<<"Cell boundary springs will be visible"<<endl;
		}else{
			cout<<"Cell boundary springs will not be visible"<<endl;
		}
		break;
	case '2':
		showPBSprings = !showPBSprings;
		if(showPBSprings){
			cout<<"Plane boundary springs will be visible"<<endl;
		}else{
			cout<<"Plane boundary springs will not be visible"<<endl;
		}
		break;
	case '3':
		showCells = !showCells;
		if(showCells){
			cout<<"Cells will be visible"<<endl;
		}else{
			cout<<"Cells will not be visible"<<endl;
		}
		break;
	case '\b':	// backspace
		init();		// reset
		break;
	case 'q':
	case 27:		// escape button
		exit(0);
		break;
	}
	glutPostRedisplay();
}

void refreshMousePos(int x, int y)
{
	double modelMatrix[16];
	double projMatrix[16];
	int    viewport[4];
	double posz; // not needed here in 2D
	
 	glGetDoublev ( GL_MODELVIEW_MATRIX,  modelMatrix );
	glGetDoublev ( GL_PROJECTION_MATRIX, projMatrix  );
	glGetIntegerv( GL_VIEWPORT,			 viewport	 );           
	gluUnProject ( x, viewport[3]-y, 0, modelMatrix, projMatrix, viewport, 
				   &mousePos.x, &mousePos.y, &posz   );
}

void motion(int x, int y)
{	// is called when the mouse is moved while button is pressed
	refreshMousePos(x,y);
	glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
	mouse_state		= state;
	mouse_button	= button;
	refreshMousePos(x,y);
	if((button == GLUT_LEFT_BUTTON)&&(state == GLUT_DOWN)){
		begin	= container.leftDown(x,y);
		end		= begin;
		if((mode==CELL)||(mode==EAST)||(mode==WEST)){	floating = false;
		}else{			floating = true; }
		switch(mode){
		case EAST:
			container.setEastCorner(begin);
			break;
		case WEST:
			container.setWestCorner(begin);
			break;
		}
		glutPostRedisplay();		
	}
	if((button == GLUT_LEFT_BUTTON)&&(state == GLUT_UP)){
		end		= container.leftDown(x,y);
		mousePos= container.getCornerPos(end);
		switch(mode){
		case CBSPRING:
			container.setSpring(begin, end);
			break;
		case PBSPRING:
			container.setSecSpring(begin, end);
			break;
		case CELL:
			motion(x,y); // update the mousePos
			container.setCell(mousePos);
			break;
		case NORTH_AND_SOUTH:
			container.setNorthAndSouth(begin, end);
			break;
		case CBTRANS:
			container.setStressSpring(begin, end);
			break;
		case PBTRANS:
			container.setSecStressSpring(begin, end);
			break;		
		}
		floating = false;
		glutPostRedisplay();
	}
}

int main(int argc, char** argv)
{
  glutInit(&argc, argv);

  glutInitWindowSize(windowWidth, windowHeight);
  glutInitWindowPosition(500,10);
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("Pick");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
//glutPassiveMotionFunc(passive);
  init();
  glutMainLoop();
  return 0;
}

/*
void passive(int x, int y) 
{

}
*/

