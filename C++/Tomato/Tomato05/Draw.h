#ifndef DRAW_H													
#define DRAW_H
//////////////////////////////////////////////////////////////////

// maximal dimension is 10, most beautifull are 9, 7, 5
#define POINT_SIZE 3
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
		if(radius>1) {
			glutWireSphere(radius, GLint(60*radius), 2);
		} else {
			glutWireSphere(radius, 20, 2);
		}
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

void drawPoint(float x,float y, float z) 
{
	glBegin(GL_POINTS); glVertex3f(x,y,z); glEnd();
}

void drawPoint(CTriple P) 
{
	drawPoint(P.x,P.y,P.z);
}

void drawDisc(CTriple P, double radius) 
{
	glPushMatrix();
		glTranslatef(P.x, P.y, P.z);
		if(radius>1) {
			glutSolidSphere(radius, GLint(60*radius), 2);
		} else {
			glutSolidSphere(radius, 20, 2);
		}
	glPopMatrix();
}

void drawLine(double x1,double y1, double z1, double x2, double y2, double z2)
{
	glBegin(GL_LINES); 
		glVertex3f(x1, y1, z1); 
		glVertex3f(x2, y2, z2); 
	glEnd();
}

void drawLine(CTriple P1, CTriple P2)
{
	drawLine(P1.x, P1.y, P1.z, P2.x, P2.y, P2.z);
}

void drawPlane(CTriple A, CTriple B, CTriple P)
{	// A x B wijst naar buiten
	// P is the beginposition where A and B originate
	// then A and B span up the plane from this point out
	glBegin(GL_QUADS);
		glVertex3f(P.x, P.y, P.z);
		glVertex3f(P.x+A.x, P.y+A.y, P.z+A.z);
		glVertex3f(P.x+A.x+B.x, P.y+A.y+B.y, P.z+A.z+B.z);
		glVertex3f(P.x+B.x, P.y+B.y, P.z+B.z);
	glEnd();
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of Draw.h												//
//////////////////////////////////////////////////////////////////