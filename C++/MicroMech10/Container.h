//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "Corner.h"

extern double gamma;		// from main.cpp
extern double floorHeight;	// from main.cpp

class CContainer
/*
 * Contains an array of the corners which are linked to one another
 * each corner has its own dynamics which is determined by spring forces
 * connecting it to its two neighbours and a pressure force
 */
{	
public:	
	CContainer					( int newNrCorners, double newRadiusCircle	);
	~CContainer					(											){};
	void	initialize			( double newCornerMass, double newPressure	);
	void	draw				(											);
	void	move				( double dt									);
	void	calculateForces		( double k, double c						);
	void	freeze				(											);
	double	calculatePressure	(											);
	double	calculateSurface	(											);
private:
	CCorner	*cornerA;
	int		nrCorners;
	double	radiusCircle;
	double	ratioPxV;
	double	pressure;
};

// constructor, initialize() must be called as well!!!
// ---------------------------------------------------
CContainer::CContainer(int newNrCorners, double newRadiusCircle)
{
	nrCorners = newNrCorners;
	radiusCircle = newRadiusCircle;
	cornerA = new CCorner[nrCorners];
}

// calculate the surface, supposing each next triangle contributes a positive part
// -------------------------------------------------------------------------------
double CContainer::calculateSurface()
{
	double	surface	= 0;	// the total surface	
	double	a,b,c;			// the length of the opposing sides
	double	s;				// half the perimeter
	CPair	A, B, C;
	
	A = cornerA[0].getPos(); 
	for(int i=1; i<nrCorners-1; ++i)
	{
		B = cornerA[i].getPos();
		C = cornerA[i+1].getPos();
		a = Distance(B,C);
		b = Distance(A,C);
		c = Distance(A,B);
		s = 0.5*(a+b+c);
		// the formula of Heron
		surface += sqrt(s*(s-a)*(s-b)*(s-c));
	}
	return surface;
}

// to freeze the current configuration
// -----------------------------------
void CContainer::freeze()
{
	cornerA[0].freeze(Distance(cornerA[nrCorners-1].getPos(),cornerA[0].getPos()));
	for(int i=1; i<nrCorners; ++i)
	{
		cornerA[i].freeze(Distance(cornerA[i-1].getPos(),cornerA[i].getPos()));
	}
}

// set the beginconditions
// -----------------------
void CContainer::initialize(double newCornerMass, double newPressure)
/*
 * Calculate beginpositions of an n-polygon of walls and corners
 * around the origin
 *
 *             corner[1]
 *                   
 *                 ^
 *               /   \
 *    wall[1]  /       \   wall[0]
 *           /           \  
 *         /               \  
 *        -------------------
 * corner[2]     wall[2]     corner[0]
 */
{
	double	angle;
	CPair	P;

	for(int i=0; i<nrCorners; ++i)
	{
		angle	= 2*PI*i/(double)(nrCorners);
		P.x		= radiusCircle*cos(angle);
		P.y		= radiusCircle*sin(angle);
		cornerA[i].initialize(P,newCornerMass);
	}
	this->freeze();
	pressure = newPressure;
	ratioPxV = pressure*this->calculateSurface();
}

// calculate a measure for the pressure
// ------------------------------------
double CContainer::calculatePressure()
{
	// suppose the surface is proportional to the volume
	return ratioPxV/this->calculateSurface();
}

// draw the current configuration
// ------------------------------
void CContainer::draw()
{	
	// the floor in blue at floorHeight
	glColor3f(0,0,1);
	drawLine(-2,floorHeight,2,floorHeight);

	// the springs connecting the corners in red
	glColor3f(1,0,0);
	for(int i=0; i<nrCorners-1; ++i)
	{
		drawLine(cornerA[i].getPos(),cornerA[i+1].getPos());
	}
	drawLine(cornerA[nrCorners-1].getPos(),cornerA[0].getPos());

	// the corners in purple
	glColor3f(1,0,1);
	for(i=0; i<nrCorners; ++i)
	{
		cornerA[i].draw();
	}
}

// move the corners by integrating the equation of motion
// ------------------------------------------------------
void CContainer::move(double dt)
/*
 * The corners updated their accelerations in the calculateForces step
 * now their positions will be updated
 */
{
	for(int i=0; i<nrCorners; ++i)
	{
		cornerA[i].move(dt);
	}
}

// calculate all the forces acting on the corners
// ----------------------------------------------
void CContainer::calculateForces(double k, double c)
/*
 * Calculate the force and modify the accelerations of the corners
 * N points from the begincorner to the endcorner (is CCW) and the
 * resulting F acts on the respective corners
 */
{
	// A // the corners are connected by springs
	/////// ------------------------------------
	CPair N;
	CPair F;
	CPair relVel;
	double distance; 
	for(int i=0; i<nrCorners-1; ++i)
	{
		N = cornerA[i+1].getPos() - cornerA[i].getPos();
		distance = VectorLength(N);
		N = N/distance; // unit length now
		F = CPair(0,0);
		relVel = cornerA[i+1].getVel() - cornerA[i].getVel();
		F = cornerA[i].calculateF(k, c, distance, N, relVel);

		// begin draw F
			CPair midPoint = cornerA[i].getPos();
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint+0.005*F);
			midPoint = cornerA[i+1].getPos();
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint-0.005*F);
		// end draw F

		cornerA[i].addF(F);
		cornerA[i+1].addF(-F);
	}
	N = cornerA[0].getPos() - cornerA[nrCorners-1].getPos();
	distance = 0;
	distance = VectorLength(N);
	N = N/distance; // unit length now
	relVel = cornerA[0].getVel() - cornerA[nrCorners-1].getVel();
	F = cornerA[nrCorners-1].calculateF(k, c, distance, N, relVel);

	// begin draw F
		CPair midPoint = cornerA[nrCorners-1].getPos();
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint+0.005*F);
		midPoint = cornerA[0].getPos();
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint-0.005*F);
	// end draw F

	cornerA[nrCorners-1].addF(F);
	cornerA[0].addF(-F);

	// B // contact with a floor on y = floorHeight
	/////// ---------------------------------------
	for(i=0; i<nrCorners; ++i)
	{
		distance = cornerA[i].getPos().y - floorHeight;
		if(distance<0)
		{
			N.x = 0; N.y = 1;
			F = cornerA[i].calculateF(k, c, distance, N, cornerA[i].getVel());

			cornerA[i].addF(-F);
		}
	}

	// C // forces exerted by the internal pressure
	/////// ---------------------------------------
	double pressure = this->calculatePressure();
	for(i=0; i<nrCorners-1; ++i)
	{
		N = cornerA[i+1].getPos() - cornerA[i].getPos();
		distance = VectorLength(N);
		N = N/distance; // unit length now
		F = CPair(N.y,-N.x); // rotation 90 deg to the right->pointing outward
		F = pow(pressure,gamma)*distance*F;

		// begin draw F
			CPair midPoint = cornerA[i].getPos()+N*0.5*distance;
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint+0.005*F);
		// end draw F

		cornerA[i].addF(F);
		cornerA[i+1].addF(F);
	}
	N = cornerA[0].getPos() - cornerA[nrCorners-1].getPos();
	distance = VectorLength(N);
	N = N/distance;
	F = CPair(N.y,-N.x); // rotation 90 deg to the right->pointing outward
	F = pow(pressure,gamma)*distance*F;

	// begin draw F
		midPoint = cornerA[nrCorners-1].getPos()+N*0.5*distance;
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint+0.005*F);
	// end draw F
	
	cornerA[nrCorners-1].addF(F);
	cornerA[0].addF(F);
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////