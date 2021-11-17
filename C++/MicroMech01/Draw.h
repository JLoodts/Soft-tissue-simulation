#ifndef DRAW_H													
#define DRAW_H
//////////////////////////////////////////////////////////////////
#define PI 3.1415926535898

// maximal dimension is 10, most beautifull are 9, 7, 5
#define POINT_SIZE 9
#define LINE_WIDTH 3.0

void drawPoint(float x,float y) 
{
	glBegin(GL_POINTS); glVertex2f(x,y); glEnd();
}

void drawDiscs(float **disc, int d, float rd) 
{
	for(int i=0; i<d; i++)
	{
		glPushMatrix();
		glTranslatef(disc[i][0], disc[i][1], 0);
		glutWireSphere(rd, GLint(60*rd), 2);
		glPopMatrix();
	}
}

void drawWalls(float **wall, int w) 
{
	glBegin(GL_LINES); 
	for(int i=0; i<w; i++)
	{
		glVertex2f(wall[i][0], wall[i][1]); 
		glVertex2f(wall[i][2], wall[i][3]); 
	}
	glEnd();
}


void drawSprings(float **wall, int w) 
{
	glBegin(GL_LINES); 
	for(int i=0; i<w-1; i++)
	{
		glVertex2f(wall[i][2], wall[i][3]); 
		glVertex2f(wall[i+1][0], wall[i+1][1]); 
	}
	glVertex2f(wall[w-1][2], wall[w-1][3]); 
	glVertex2f(wall[0][0], wall[0][1]); 
	glEnd();
}

//////////////////////////////////////////////////////////////////
#endif