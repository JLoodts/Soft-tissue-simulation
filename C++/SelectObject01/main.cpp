/* 
 *  Demonstration of picking and rendering luminous objects.  Drag the
 *  middle mouse button to spin the object.  Move the mouse over the
 *  bulbs to light them.
 * 
 *  author: Nate Robins
 *  email: ndr@pobox.com
 *  www: http://www.pobox.com/~ndr 
 */


/* includes */
#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include <GL/glut.h>


/* defines */
//#define SELECT_BUFFER 32


/* globals */
GLuint signal_list = 0;			/* display list for traffic signal */

const GLsizei SELECT_BUFFER = 32;
GLuint select_buffer[SELECT_BUFFER];	/* selection buffer */
GLuint selected;	// 0 : no particle is selected else the number of the light
bool selectOne = false;
GLuint picked = 0;			/* current light that is picked */			
double previousMousePosX = 0;
double previousMousePosY = 0;

GLint mouse_state = -1;
GLint mouse_button = -1;



GLuint windowWidth = 512;			/* in pixels */
GLuint windowHeight = windowWidth;

GLUquadricObj* quadric;

double** pos;


/* functions */



void
init(void)
{

  /* lighting */
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  glEnable(GL_CULL_FACE);

  glClearColor(0.3, 0.3, 1.0, 0.0);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  /* set the selection buffer size */
  glSelectBuffer(SELECT_BUFFER, select_buffer);

  pos	= new double*[3];
  for(int i=0; i<3; ++i){ pos[i] = new double[3];} 
  pos[0][0] = 0.0; // red light
  pos[0][1] = 1.75;
  pos[0][2] = -8.0;
  pos[1][0] = 0.0;	// yellow light
  pos[1][1] = 0.0;
  pos[1][2] = -8.2;
  pos[2][0] = 0.0;	// green light
  pos[2][1] = -1.75;
  pos[2][2] = -8.1;
}

void
reshape(int width, int height)
{
  GLfloat light_position[4] = { 1.0, 1.0, 1.0, 0.0 };

  glViewport(0, 0, width, height);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (GLfloat)height / (GLfloat)width, 1.0, 128.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}

void
drawObjects(void)
{
	glPushMatrix();
	
//	glTranslatef(0.0, 0.0, -8.0);
//	glRotatef(15.0, 0.0, 1.0, 0.0);
	
	

//	glEnable(GL_LIGHT1);

	/* draw the lights */
	glPushMatrix();
	glTranslatef(pos[0][0], pos[0][1], pos[0][2]);
	if (selected == 1) { // red light
		glDisable(GL_LIGHTING);
		glColor3f(1.0, 0.0, 0.0);
		glLoadName(1);
		glutSolidSphere(0.75, 32, 32);
		glEnable(GL_LIGHTING);
	} else {
		glColor3f(0.3, 0.0, 0.0);
		glLoadName(1);
		glutSolidSphere(0.75, 32, 32);
	}
	glPopMatrix();
	glPushMatrix();
	glTranslatef(pos[1][0], pos[1][1], pos[1][2]);
	if (selected == 2) { // yellow light
		glDisable(GL_LIGHTING);
		glColor3f(1.0, 1.0, 0.0);
		glLoadName(2);
		glutSolidSphere(0.75, 32, 32);
		glEnable(GL_LIGHTING);
	} else {
		glColor3f(0.3, 0.3, 0.0);
		glLoadName(2);
		glutSolidSphere(0.75, 32, 32);
	}
	glPopMatrix();
	glPushMatrix();
	glTranslatef(pos[2][0], pos[2][1], pos[2][2]);
	if (selected == 3) { // green light
		glDisable(GL_LIGHTING);
		glColor3f(0.0, 1.0, 0.0);
		glLoadName(3);
		glutSolidSphere(0.75, 32, 32);
		glEnable(GL_LIGHTING);
	} else {
		glColor3f(0.0, 0.3, 0.0);
		glLoadName(3);
		glutSolidSphere(0.75, 32, 32);
	}
	glPopMatrix();

	glPopMatrix();
}

void
display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  drawObjects();

  glutSwapBuffers();
}

void
keyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'h':
    printf("help\n\n");
    printf("f            -  Filled\n");
    printf("w            -  Wireframe\n");
    printf("backspace    -  Reset\n");
    printf("escape or q  -  Quit\n\n");
    break;

  case 'f':
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    break;

  case 'w':
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    break;

  case 's':
    selectOne = true;
	cout<<"Left click the light you want to select."<<endl;
	cout.flush();
    break;

  case '\b':
    init();
    break;

  case 'q':
  case 27:
    exit(0);
    break;
  }

  glutPostRedisplay();
}

GLuint
pick(int x, int y)
{
  GLuint    i, hits, num_names, picked;
  GLuint*   p;
  GLboolean save;
  GLint     viewport[4];
  GLuint    depth = (GLuint)-1;

  /* get the current viewport parameters */
  glGetIntegerv(GL_VIEWPORT, viewport);

  /* set the render mode to selection */
  glRenderMode(GL_SELECT);
  
  glInitNames(); // clears the name stack
  glPushName(0); // puts 0 on top of the name stack
  //GLuint d[1]; glGetIntegerv(GL_NAME_STACK_DEPTH, d); //-> 3.435.97.3836 names possible

  /* setup a picking matrix and render into selection buffer */
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();

  glLoadIdentity();
  gluPickMatrix(x, viewport[3] - y, 1.0, 1.0, viewport); 
  // must be called BEFORE gluPerspective() or glOrtho()
  gluPerspective(60.0, (GLfloat)viewport[3]/(GLfloat)viewport[2], 1.0, 128.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  drawObjects();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  hits = glRenderMode(GL_RENDER);

  p = select_buffer;
  picked = 0;
  for (i = 0; i < hits; i++) {
    save = GL_FALSE;
    num_names = *p;			/* number of names in this hit */
    p++;

    if (*p <= depth) {			/* check the 1st depth value */
      depth = *p;
      save = GL_TRUE;
    }
    p++;
    if (*p <= depth) {			/* check the 2nd depth value */
      depth = *p;
      save = GL_TRUE;
    }
    p++;

    if (save)
      picked = *p;

    p += num_names;			/* skip over the rest of the names */
  }

  return picked;
}

double calculatePointDepth(const double &x, const double &y, const double &z)
 {
  // this figures our where inbetween the near and far clipping planes a point is.
  // use the result of this with the Pick function so that the mouse is the same
  // depth into the screen as that point.

  double modelMatrix[16];
  double projMatrix[16];
  double junkx, junky, depth;
  int    viewport[4];
                
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  gluProject(x, y, z, modelMatrix, projMatrix, viewport, &junkx, &junky, &depth);
  
  return depth;
 }

void
mouse(int button, int state, int x, int y)
{
	mouse_state = state;
	mouse_button = button;
	previousMousePosX = x;
	previousMousePosY = y;
	if((selectOne)&&(button == GLUT_LEFT_BUTTON)&&(state == GLUT_DOWN)){	
		selected = picked;
		selectOne = false;
	}
}

void // is called when the mouse is moved while button is pressed
motion(int x, int y)
{
	if(selected>0) {
		double modelMatrix[16];
		double projMatrix[16];
		int    viewport[4];
		double posx, posy, posz;
        double pointDepth = calculatePointDepth(pos[selected-1][0],pos[selected-1][1],pos[selected-1][2]);     
		glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
		glGetIntegerv(GL_VIEWPORT, viewport);           
		gluUnProject(x, viewport[3]-y, pointDepth, modelMatrix, projMatrix, viewport, &posx, &posy, &posz);
		pos[selected-1][0] = posx;
		pos[selected-1][1] = posy;
		pos[selected-1][2] = posz; 
		glutPostRedisplay();
	}
}

void
passive(int x, int y) 
{
  picked = pick(x, y);
  glutPostRedisplay();
}

int
main(int argc, char** argv)
{
  glutInit(&argc, argv);

  glutInitWindowSize(windowWidth, windowHeight);
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("Pick");
  
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutPassiveMotionFunc(passive);
  
  init();
  
  glutMainLoop();
  return 0;
}
