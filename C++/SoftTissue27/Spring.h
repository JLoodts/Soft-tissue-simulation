//////////////////////////////////////////////////////////////////
#ifndef SPRING_H												//
#define SPRING_H												//
// begin of Spring.h											//
//////////////////////////////////////////////////////////////////

#include "Corner.h"

extern double	scaleFactorForForces;	// from main.cpp
extern bool		showForces;		// from main.cpp

/**
 * Spring class	
 * models a piece of cell wall, its features are:
 *  - a linear spring connecting two corners 
 *	- restLength keeps track of the length of the spring in rest
 */
class Spring
{	
public:	
	Spring(){}
	~Spring(){}
	void	initialize( int newNr, double newK, double newC, 
					    Corner* newBegin, Corner* newEnd );
	void	initialize( int newNr, double newK, double newC, 
					    Corner* newBegin, Corner* newEnd, double newrestLength );
	double	getLoad(){ return load; }
	void	draw();
	void	calculateForce();
	int		getNr(){ return nr; }
	void	save( ofstream outFile ){	outFile<<nr<<" "<<pBegin->getNr()<<" "<<pEnd->getNr()<<endl; }
private:
	void	showForce( CTriple F );
	void	addForce( CTriple F );
	int		nr;
	Corner*	pBegin;			// position of the begin point
	Corner*	pEnd;			// position of the end point
	double	restLength;			// the length of the spring in rest
	double	load;				// the current load (N) on this spring
	double	k;
	double	c;
};

void Spring::initialize(int newNr, double newK, double newC, Corner* newBegin, Corner* newEnd)
{
	this->initialize(newNr,newK,newC,newBegin,newEnd,Distance(newBegin->getPos(),newEnd->getPos()));
}

void Spring::initialize(int newNr, double newK, double newC, Corner* newBegin, Corner* newEnd, double newRestLength)
{
	nr			= newNr;
	k			= newK;
	c			= newC;
	pBegin		= newBegin;
	pEnd		= newEnd;
	restLength	= newRestLength;
	load		= 0; // this may not be true since the spring might be deformed already
}

void Spring::draw()
{
//	double color = (load)/(loadMax);
//	glColor3f(color,0,1-color); // red = max, blue = min
//	glColor3f(1,0,0);
	drawLine(pBegin->getPos(),pEnd->getPos());
}

void Spring::showForce(CTriple F)
{
	CTriple midPoint = pBegin->getPos();
	glColor3f(1,0,0);
	drawLine(midPoint ,midPoint+scaleFactorForForces*F);
	midPoint = pEnd->getPos();
	glColor3f(1,0,0);
	drawLine(midPoint ,midPoint-scaleFactorForForces*F);
}

void Spring::addForce(CTriple F)
{
	pBegin->addF(F);
	pEnd->addF(-F);
}

void Spring::calculateForce()
{
	CTriple N, relVel, F, midPoint;
	double dx, distance;

	relVel = pEnd->getVel() - pBegin->getVel();
	N = pEnd->getPos() - pBegin->getPos();
	distance = VectorLength(N);
	dx = restLength - distance;
	if(dx<0){
		N = N/distance; // unit length now
		double	dv = relVel*N;
		F = -(k*dx-c*dv)*N;
	} else {
		F = ZEROTRIPLE;
	}
	
	load = fabs(F.x);
//	if(load>loadMax){ loadMax = load; }

	if(showForces){	showForce(F); }
	
	addForce(F);
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Spring.h												//
//////////////////////////////////////////////////////////////////