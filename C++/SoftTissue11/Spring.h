//////////////////////////////////////////////////////////////////
#ifndef SPRING_H												//
#define SPRING_H												//
// begin of Spring.h											//
//////////////////////////////////////////////////////////////////

#include "Corner.h"

extern double loadMax;		// from main.cpp
extern double scaleFactor;	// from main.cpp
extern bool	showForces;		// from main.cpp

class CSpring
/**
 *	models a corner of the polygon its features are:
 *  - an linear spring connecting two corners: distanceToNext keeps
 *    track of the length of the spring in rest
 *	- they are the only things that move
 */
{	
public:	
	CSpring(){}
	~CSpring(){}
	void	initialize(int newNr, CCorner* newBegin, CCorner* newEnd);
	void	initialize(int newNr, CCorner* newBegin, CCorner* newEnd, double newRestLength);
	int		getNr(){return nr;}
	double	getLoad(){return load;}
	double	getRestLength(){return restLength;}
	void	setRestLength(double newRestLength){restLength = newRestLength;}
	void	draw();
	void	calculateFspring(double k, double c);
private:
	void	showForce(CPair F);
	void	addForce(CPair F);
	int		nr;
	CCorner*	pBegin;			// position of the begin point
	CCorner*	pEnd;			// position of the end point
	double	restLength;			// the length of the spring in rest
	double	load;				// the current load (N) on this spring
};

void CSpring::initialize(int newNr, CCorner* newBegin, CCorner* newEnd)
{
	this->initialize(newNr,newBegin,newEnd,Distance(newBegin->getPos(),newEnd->getPos()));
}

void CSpring::initialize(int newNr, CCorner* newBegin, CCorner* newEnd, double newRestLength)
{
	nr = newNr;
	pBegin = newBegin;
	pEnd = newEnd;
	restLength = newRestLength;
	load = 0;
}

void CSpring::draw()
{
	double color = (load)/(loadMax);
	glColor3f(color,0,1-color); // red = max, blue = min
	drawLine(pBegin->getPos(),pEnd->getPos());
}

void CSpring::showForce(CPair F)
{
	CPair midPoint = pBegin->getPos();
	glColor3f(0,1,0);
	drawLine(midPoint ,midPoint+scaleFactor*F);
	midPoint = pEnd->getPos();
	glColor3f(0,1,0);
	drawLine(midPoint ,midPoint-scaleFactor*F);
}

void CSpring::addForce(CPair F)
{
	pBegin->addF(F);
	pEnd->addF(-F);
}

void CSpring::calculateFspring(double k, double c)
{
	CPair N, relVel, F, midPoint;
	double dx, distance;

	relVel = pEnd->getVel() - pBegin->getVel();
	N = pEnd->getPos() - pBegin->getPos();
	distance = VectorLength(N);
	dx = restLength - distance;
	if(dx<0){
		N = N/distance; // unit length now
		F = CPair(0,0);
		F = pBegin->calculateF(k, c, dx, N, relVel);
	} else {
		F = ZEROPAIR;
	}
	
	load = VectorLength(F);
	if(load>loadMax){
		loadMax = load;
	}

	if(showForces){
		showForce(F);
	}
	
	addForce(F);
}



//////////////////////////////////////////////////////////////////
#endif															//
// end of Spring.h												//
//////////////////////////////////////////////////////////////////