#ifndef DRAW_H													
#define DRAW_H
//////////////////////////////////////////////////////////////////

// maximal dimension is 10, most beautifull are 9, 7, 5
#define POINT_SIZE 1
#define LINE_WIDTH 1

void drawPoint(float x,float y) 
{
	glBegin(GL_POINTS); glVertex2f(x,y); glEnd();
}

void drawPoint(CPair P) 
{
	drawPoint(P.x,P.y);
}

void drawDisc(CPair P, double radius) 
{
	glPushMatrix();
		glTranslatef(P.x, P.y, 0);
		glutWireSphere(radius, GLint(60*radius), 2);
	glPopMatrix();
}

void drawLine(double x1,double y1, double x2, double y2)
{
	glBegin(GL_LINES); 
		glVertex2f(x1, y1); 
		glVertex2f(x2, y2); 
	glEnd();
}

void drawLine(CPair P1, CPair P2)
{
	drawLine(P1.x, P1.y, P2.x, P2.y);
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of Draw.h												//
//////////////////////////////////////////////////////////////////