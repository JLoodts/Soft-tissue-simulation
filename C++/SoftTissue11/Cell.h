//////////////////////////////////////////////////////////////////
#ifndef CELL_H													//
#define CELL_H													//
// begin of Cell.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "InputFile.h"
#include "Spring.h"

extern double kCorner;		// from main.cpp
extern bool	showForces;		// from main.cpp
extern double scaleFactor;	// from main.cpp
extern double kPressure;	// from main.cpp
extern int	  forceMode;	// from main.cpp

#include <fstream.h>	// for file handling

class CCell
/*
 *	models a cell which is comprised of:
 *  - an array of pointers to the cornerpoints
 *	- a measure for the internal pressure based on its surface/volume
 */
{	
public:	
	CCell(){}
	void initialize(int newNr, int newNrCorners, int* cornerNrs, CCorner* cornerA, double newPressure	);
	void setInitialConditions(double newRatioPxV, double newRestPressure){ratioPxV = newRatioPxV; restPressure = newRestPressure;}
	~CCell(){}
	void	draw		(								);
	void	drawFull();
	double	calculateSurface	(											);
	double	calculateDeltaPressure	(											);
	int		getNrCorners(){return nrCorners;}
	int		getCornerNr(int i){return pCornerA[i]->getNr();}
	void	calculateF	( double k,  double c);
//	double	getLoadInWall(int i){return springA[i].getLoad();}
	void	drawWall(int i);
	int		getNr(){return nr;}
	CSpring* getWallP(int i){return &springA[i];}
	void	save(ofstream outFile);
	void	setRestLength(int i, double length){springA[i].setRestLength(length);}
private:
	int		nr;
	CCorner** pCornerA;
	int nrCorners;
	//double pressure;
	double ratioPxV;
	double restPressure;
	CSpring*	springA;
};

// set the initial values
// ----------------------
void CCell::initialize(int newNr, int newNrCorners, int* cornerNrs, CCorner* cornerA, double newPressure)
{
	nr = newNr;
	nrCorners = newNrCorners;
	pCornerA = new CCorner*[nrCorners];
	
	for(int i=0; i<nrCorners; ++i)
	{
		pCornerA[i] = &cornerA[cornerNrs[i]-1]; 
	}
	restPressure = newPressure;
	ratioPxV = newPressure*this->calculateSurface();
	
	// set the restlength of the longitudinal springs in the cell wall
	springA = new CSpring[nrCorners];
	for(i=0; i<nrCorners-1; ++i)
	{
		springA[i].initialize(i,pCornerA[i],pCornerA[i+1]);
	}
	springA[nrCorners-1].initialize(nrCorners-1,pCornerA[nrCorners-1],pCornerA[0]);
		
}

// calculate a measure for the pressure
// ------------------------------------
double CCell::calculateDeltaPressure()
{
	// suppose the surface is proportional to the volume
	return ratioPxV/this->calculateSurface() - restPressure;
}

void CCell::drawWall(int i)
{
	springA[i].draw();
}


// calculate the force due to the spring with the next corner (in CCW direction)
// -----------------------------------------------------------------------------
void CCell::calculateF(double k, double c)
{
	int i;
	double distance;
	CPair N, relVel, F, midPoint;
	// mode = 1: everything
	//		= 2: only A

	
	// A // forces exerted by the connecting springs
	/////// ----------------------------------------
	for(i=0; i<nrCorners-1; ++i)
	{
		springA[i].calculateFspring(k,c);
	}
	springA[i].calculateFspring(k,c);	

	if(forceMode == 1) {
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

		pCornerA[i]->addFtoAccPressure(F);
		pCornerA[i+1]->addFtoAccPressure(F);
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
	
	pCornerA[nrCorners-1]->addFtoAccPressure(F);
	pCornerA[0]->addFtoAccPressure(F);
		

	} // end if(mode == 1) 

	// C // forces exerted by the big springs connecting corners with 1 corner in between
	/////// -----------------------------------------
/*	relVel = pCornerA[1]->getVel() - pCornerA[nrCorners-1]->getVel();
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
*/
}



void CCell::draw() 
{ 
	int i;
	// the springs connecting the corners in red
	//glColor3f(1,0,0);
	//double color;
/**/	for(i=0; i<nrCorners-1; ++i)
	{
		//color = 0.5*(pCornerA[i]->getLoad()+pCornerA[i+1]->getLoad())/(loadMax);
		//glColor3f(color,0,1-color); // red = max, blue = min
		//drawLine(pCornerA[i]->getPos(),pCornerA[i+1]->getPos());
		drawWall(i);
	}
	//color = 0.5*(pCornerA[nrCorners-1]->getLoad()+pCornerA[0]->getLoad())/(loadMax);
	//glColor3f(color,0,1-color); // red = max, blue = min
	//drawLine(pCornerA[nrCorners-1]->getPos(),pCornerA[0]->getPos());
	drawWall(nrCorners-1);
	
	// the corners in purple
	//glColor3f(1,0,1);
/*	for(i=0; i<nrCorners; ++i)
	{
		color = (pCornerA[i]->getLoad())/(loadMax);
		glColor3f(color,0,1-color); // red = max, blue = min
		pCornerA[i]->draw();
	}
*/
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
/*	double	surface	= 0;	// the total surface	
	double	a,b,c;			// the length of the opposing sides
	double	s;				// half the perimeter
	CPair	A, B, C;
	int i;

	A = pCornerA[0]->getPos();
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
*/
	double	surface	= 0;	// the total surface
	int i;
	CPair	A, B, C;
	double IvectorProdI;
	A = pCornerA[0]->getPos();

	for(i=0; i<nrCorners-1; ++i)
	{
		B = pCornerA[i]->getPos()-A;
		C = pCornerA[i+1]->getPos()-A;
		IvectorProdI = B.x*C.y-B.y*C.x;

		surface += 0.5*IvectorProdI;
	}

	

	return -surface;
}

void CCell::save(ofstream outFile)
{
	outFile<<"cellnr "<<nr<<endl;
	outFile<<"ratiopxv "<<ratioPxV<<endl;
	outFile<<"restpressure "<<restPressure<<endl;
	outFile<<"restlength ";
	for(int i=0; i<nrCorners; ++i)
	{
		outFile<<springA[i].getRestLength()<<" ";
	}
	outFile<<endl;
}


//////////////////////////////////////////////////////////////////
#endif															//
// end of Cell.h												//
//////////////////////////////////////////////////////////////////