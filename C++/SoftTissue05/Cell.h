//////////////////////////////////////////////////////////////////
#ifndef CELL_H													//
#define CELL_H													//
// begin of Cell.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
extern double kCorner;		// from main.cpp
extern bool	showForces;		// from main.cpp
extern double loadMax;		// from main.cpp
extern double scaleFactor;	// from main.cpp
extern double kPressure;	// from main.cpp

class CCell
/*
 *	models a cell which is comprised of:
 *  - an array of pointers to the cornerpoints
 *	- a measure for the internal pressure based on its surface/volume
 */
{	
public:	
	CCell(){}
	void initialize(int newNrCorners, int* cornerNrs, CCorner* cornerA, double newPressure	);
	~CCell(){}
	void	draw		(								);
	void	drawFull();
	double	calculateSurface	(											);
	double	calculateDeltaPressure	(											);

/*	CPair	getPos		(								) { return pos;				}
	CPair	getVel		(								) { return vel;				}
*/	void	calculateF	( double k,  double c, int mode);
	//	void	move		( double dt						);
	
	void	freeze		(); 
private:
	CCorner** pCornerA;
	int nrCorners;
	//double pressure;
	double ratioPxV;
	double restPressure;
	double	*restLengthA;	// distance to the neighbouring corners: to calculate the springforce
	double	*restLength1CornerA;	// springs connecting cornerpoints with 1 corner in between
	//double	*restLength2CornersA; // springs connecting cornerpoints with 2 corners in between
};

// set the initial values
// ----------------------
void CCell::initialize(int newNrCorners, int* cornerNrs, CCorner* cornerA, double newPressure)
{
	nrCorners = newNrCorners;
	pCornerA = new CCorner*[nrCorners];
	for(int i=0; i<nrCorners; ++i)
	{
		pCornerA[i] = &cornerA[cornerNrs[i]-1]; 
	}
	restPressure = newPressure;
	ratioPxV = newPressure*this->calculateSurface();
	
	// set the restlength of the longitudinal springs in the cell wall
	restLengthA = new double[nrCorners];
	for(i=0; i<nrCorners-1; ++i)
	{
		restLengthA[i] = Distance(pCornerA[i]->getPos(),pCornerA[i+1]->getPos());
	}
	restLengthA[nrCorners-1] = Distance(pCornerA[nrCorners-1]->getPos(),pCornerA[0]->getPos());

	// set the restlength of the springs connecting the previous and the next corner
	restLength1CornerA = new double[nrCorners];
	CPair springVector = pCornerA[1]->getPos()-pCornerA[nrCorners-1]->getPos();
	restLength1CornerA[0] = VectorLength(springVector);
	for(i=1; i<nrCorners-1; ++i)
	{
		springVector = pCornerA[i+1]->getPos()-pCornerA[i-1]->getPos();
		restLength1CornerA[i] = VectorLength(springVector);
	}
	springVector = pCornerA[0]->getPos()-pCornerA[nrCorners-2]->getPos();
	restLength1CornerA[nrCorners-1] = VectorLength(springVector);

		
}

// calculate a measure for the pressure
// ------------------------------------
double CCell::calculateDeltaPressure()
{
	// suppose the surface is proportional to the volume
	return ratioPxV/this->calculateSurface() - restPressure;
}

// calculate the force due to the spring with the next corner (in CCW direction)
// -----------------------------------------------------------------------------
void CCell::calculateF(double k, double c, int mode)
{
	int i;
	double dx, distance;
	CPair N, relVel, F, midPoint;
	// mode = 1: everything
	//		= 2: only A

	
	// A // forces exerted by the connecting springs
	/////// ----------------------------------------
	for(i=0; i<nrCorners-1; ++i)
	{
		relVel = pCornerA[i+1]->getVel() - pCornerA[i]->getVel();
		N = pCornerA[i+1]->getPos() - pCornerA[i]->getPos();
		distance = VectorLength(N);
		dx = restLengthA[i] - distance;
		N = N/distance; // unit length now
		F = CPair(0,0);
		F = pCornerA[i]->calculateF(k, c, dx, N, relVel);
		if(showForces){// begin draw F
			midPoint = pCornerA[i]->getPos();
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint+scaleFactor*F);
			midPoint = pCornerA[i+1]->getPos();
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint-scaleFactor*F);
		// end draw F
		}
		pCornerA[i]->addF(F);
		pCornerA[i+1]->addF(-F);
	}
	relVel = pCornerA[0]->getVel() - pCornerA[nrCorners-1]->getVel();
	N = pCornerA[0]->getPos() - pCornerA[nrCorners-1]->getPos();
	distance = VectorLength(N);
	dx = restLengthA[nrCorners-1] - distance;
	N = N/distance; // unit length now
	F = CPair(0,0);
	F = pCornerA[nrCorners-1]->calculateF(k, c, dx, N, relVel);
	if(showForces){// begin draw F
		midPoint = pCornerA[nrCorners-1]->getPos();
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint+scaleFactor*F);
		midPoint = pCornerA[0]->getPos();
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint-scaleFactor*F);
	// end draw F
	}
	pCornerA[nrCorners-1]->addF(F);
	pCornerA[0]->addF(-F);

	
	if(mode == 1) {
	// B // forces exerted by the internal pressure
	/////// ---------------------------------------

	double dPressure = this->calculateDeltaPressure();
	for(i=0; i<nrCorners-1; ++i)
	{
		N = pCornerA[i+1]->getPos() - pCornerA[i]->getPos();
		distance = VectorLength(N);
		N = N/distance; // unit length now
		F = CPair(-N.y,N.x); // rotation 90 deg to the left->pointing outward
		F = kPressure*dPressure*distance*F;

		if(showForces){// begin draw F
			CPair midPoint = pCornerA[i]->getPos()+N*0.5*distance;
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint+scaleFactor*F);
		// end draw F
		}

		pCornerA[i]->addF(F);
		pCornerA[i+1]->addF(F);
	}
	N = pCornerA[0]->getPos() - pCornerA[nrCorners-1]->getPos();
	distance = VectorLength(N);
	N = N/distance;
	F = CPair(-N.y,N.x); // rotation 90 deg to the left->pointing outward
	F = kPressure*dPressure*distance*F;

	if(showForces){// begin draw F
		midPoint = pCornerA[nrCorners-1]->getPos()+N*0.5*distance;
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint+scaleFactor*F);
	// end draw F
	}
	
	pCornerA[nrCorners-1]->addF(F);
	pCornerA[0]->addF(F);

	} // end if(mode == 1) 

	// C // forces exerted by the big springs connecting corners with 1 corner in between
	/////// -----------------------------------------
	relVel = pCornerA[1]->getVel() - pCornerA[nrCorners-1]->getVel();
	N = pCornerA[1]->getPos()-pCornerA[nrCorners-1]->getPos();
	distance = VectorLength(N);
	dx = restLength1CornerA[0] - distance;
	N = N/distance; // unit length now
	F = pCornerA[0]->calculateF(kCorner, c, dx, N, relVel);
	if(showForces){// begin draw F
		midPoint = pCornerA[nrCorners-1]->getPos();
		glColor3f(0,1,1);
		drawLine(midPoint ,midPoint+scaleFactor*F);
		midPoint = pCornerA[1]->getPos();
		glColor3f(0,1,1);
		drawLine(midPoint ,midPoint-scaleFactor*F);
	// end draw F
	}
	pCornerA[nrCorners-1]->addF(F);
	pCornerA[1]->addF(-F);
	for(i=1; i<nrCorners-1; ++i)
	{
		relVel = pCornerA[i+1]->getVel() - pCornerA[i-1]->getVel();
		N = pCornerA[i+1]->getPos()-pCornerA[i-1]->getPos();
		distance = VectorLength(N);
		dx = restLength1CornerA[i] - distance;
		N = N/distance; // unit length now
		F = pCornerA[i]->calculateF(kCorner, c, dx, N, relVel);
		if(showForces){// begin draw F
			midPoint = pCornerA[i-1]->getPos();
			glColor3f(0,1,1);
			drawLine(midPoint ,midPoint+scaleFactor*F);
			midPoint = pCornerA[i+1]->getPos();
			glColor3f(0,1,1);
			drawLine(midPoint ,midPoint-scaleFactor*F);
		// end draw F
		}
		pCornerA[i-1]->addF(F);
		pCornerA[i+1]->addF(-F);
	}
	relVel = pCornerA[0]->getVel() - pCornerA[nrCorners-2]->getVel();
	N = pCornerA[0]->getPos()-pCornerA[nrCorners-2]->getPos();
	distance = VectorLength(N);
	dx = restLength1CornerA[nrCorners-1] - distance;
	N = N/distance; // unit length now
	F = pCornerA[nrCorners-1]->calculateF(kCorner, c, dx, N, relVel);
	if(showForces){// begin draw F
		midPoint = pCornerA[nrCorners-2]->getPos();
		glColor3f(0,1,1);
		drawLine(midPoint ,midPoint+scaleFactor*F);
		midPoint = pCornerA[0]->getPos();
		glColor3f(0,1,1);
		drawLine(midPoint ,midPoint-scaleFactor*F);
	// end draw F
	}
	pCornerA[nrCorners-2]->addF(F);
	pCornerA[0]->addF(-F);

}


void CCell::freeze()
{
	// set the restlength
	restLengthA = new double[nrCorners];
	for(int i=0; i<nrCorners-1; ++i)
	{
		restLengthA[i] = Distance(pCornerA[i]->getPos(),pCornerA[i+1]->getPos());
	}
	restLengthA[nrCorners-1] = Distance(pCornerA[nrCorners-1]->getPos(),pCornerA[0]->getPos());

	// set the restlength of the springs connecting the 3/2th point of the wall and the middle of the next one
	CPair springVector = pCornerA[1]->getPos()-pCornerA[nrCorners-1]->getPos();
	restLength1CornerA[0] = VectorLength(springVector);
	for(i=1; i<nrCorners-1; ++i)
	{
		springVector = pCornerA[i+1]->getPos()-pCornerA[i-1]->getPos();
		restLength1CornerA[i] = VectorLength(springVector);
	}
	springVector = pCornerA[0]->getPos()-pCornerA[nrCorners-2]->getPos();
	restLength1CornerA[nrCorners-1] = VectorLength(springVector);

}

void CCell::draw() 
{ 
	int i;
	// the springs connecting the corners in red
	//glColor3f(1,0,0);
	double color;
/**/	for(i=0; i<nrCorners-1; ++i)
	{
		color = 0.5*(pCornerA[i]->getLoad()+pCornerA[i+1]->getLoad())/(loadMax);
		glColor3f(color,0,1-color); // red = max, blue = min
		drawLine(pCornerA[i]->getPos(),pCornerA[i+1]->getPos());
	}
	color = 0.5*(pCornerA[nrCorners-1]->getLoad()+pCornerA[0]->getLoad())/(loadMax);
	glColor3f(color,0,1-color); // red = max, blue = min
	drawLine(pCornerA[nrCorners-1]->getPos(),pCornerA[0]->getPos());
	
	// the corners in purple
	//glColor3f(1,0,1);
	for(i=0; i<nrCorners; ++i)
	{
		color = (pCornerA[i]->getLoad())/(loadMax);
		glColor3f(color,0,1-color); // red = max, blue = min
		pCornerA[i]->draw();
	}
}

void CCell::drawFull() 
{ 
	int i;
	glColor3f(1,0,0);
	glBegin(GL_POLYGON);
		for(i=0; i<nrCorners; ++i)
		{
			glVertex2f(pCornerA[i]->getPos().x,pCornerA[i]->getPos().y);
		}
	glEnd();
}

// calculate the surface, supposing each next triangle contributes a positive part
// -------------------------------------------------------------------------------
double CCell::calculateSurface()
{
	double	surface	= 0;	// the total surface	
	double	a,b,c;			// the length of the opposing sides
	double	s;				// half the perimeter
	CPair	A, B, C;
	int i;
/*	A = ZEROPAIR;
	for(i=0; i<nrCorners-1; ++i)
	{ A = A + pCornerA[i]->getPos(); }
	A = A/nrCorners; 
	glColor3f(0,1,0);
	drawPoint(A.x,A.y);
	for(i=0; i<nrCorners-1; ++i)
	{
		B = pCornerA[i]->getPos();
		C = pCornerA[i+1]->getPos();
		a = Distance(B,C);
		b = Distance(A,C);
		c = Distance(A,B);
		s = 0.5*(a+b+c);
		// the formula of Heron
		surface += sqrt(s*(s-a)*(s-b)*(s-c));
	}*/
/**/	A = pCornerA[0]->getPos();
	for(i=0; i<nrCorners-1; ++i)
	{
		B = pCornerA[i]->getPos();
		C = pCornerA[i+1]->getPos();
		a = Distance(B,C);
		b = Distance(A,C);
		c = Distance(A,B);
		s = 0.5*(a+b+c);
		// the formula of Heron
		surface += sqrt(s*(s-a)*(s-b)*(s-c));
	}
	
	return surface;
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of Cell.h												//
//////////////////////////////////////////////////////////////////