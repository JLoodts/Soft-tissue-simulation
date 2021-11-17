#ifndef DRAW_H													
#define DRAW_H
//////////////////////////////////////////////////////////////////

// maximal dimension is 10, most beautifull are 9, 7, 5
#define POINT_SIZE 1
#define LINE_WIDTH 1.0

void DrawPoint(float x,float y) 
{
	glBegin(GL_POINTS); glVertex2f(x,y); glEnd();
}

void DrawPoint(CPair P) 
{
	DrawPoint(P.x,P.y);
}

void DrawDisc(CPair P, double radius) 
{
	glPushMatrix();
		glTranslatef(P.x, P.y, 0);
		glutSolidSphere(radius, GLint(30), 2);
	glPopMatrix();
}

void DrawLine(double x1,double y1, double x2, double y2)
{
	glBegin(GL_LINES); 
		glVertex2f(x1, y1); 
		glVertex2f(x2, y2); 
	glEnd();
}

void DrawLine(CPair P1, CPair P2)
{
	DrawLine(P1.x, P1.y, P2.x, P2.y);
}

void DrawRectangle(double LowLeftx, double LowLefty, double UpRightx, double UpRighty)
{
	glRectf(LowLeftx, LowLefty, UpRightx, UpRighty);
}

void DrawRectangle(CPair LowLeft, CPair UpRight)
{
	DrawRectangle(LowLeft.x, LowLeft.y, UpRight.x, UpRight.y);
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of Draw.h												//
//////////////////////////////////////////////////////////////////