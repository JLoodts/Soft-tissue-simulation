//////////////////////////////////////////////////////////////////
#ifndef CORNER_H												//
#define CORNER_H												//
// begin of Corner.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
extern double g;

class CCorner
/*
 *	models a corner of the polygon its features are:
 *  - an linear spring connecting two corners: distanceToNext keeps
 *    track of the length of the spring in rest
 *	- they are the only things that move
 */
{	
public:	
	CCorner(){}
	~CCorner(){}
	void	initialize	( CPair newPos, double newMass	);
	CPair	getPos		(								) { return pos;				}
	void	setPosY		( double newY					) { pos.y = newY;			}
	CPair	getVel		(								) { return vel;				}
	void	setVel		( CPair newVel					) { vel = newVel;			}
	void	setAcc		( CPair newAcc					) { acc = newAcc;			}
	double	getLoad		(								) { return load;			}
	CPair	calculateF	( double k,  double c, double distance, CPair N, CPair relVel );
	void	draw		(								) { drawPoint(pos);			}
	void	move		( double dt						);
	void	addF		( CPair &F						) { acc = acc + (F/mass); load += VectorLength(F);	}
	CPair	getF		(								) { return mass*acc;		}
private:
	CPair	pos;			// position
	CPair	vel;			// velocity
	CPair	acc;			// acceleration
	double	mass;			// mass of a corner
	double	load;			// sum of magnitude of the forces
};

// set the initial values
// ----------------------
void CCorner::initialize(CPair newPos, double newMass)
{
	pos		= newPos;
	vel		= ZEROPAIR;
	acc		= ZEROPAIR;
	mass	= newMass;
	load	= 0.0;
	// freeze must be called now in Container in order to set distanceToNext
}

// calculate the force due to the spring with the next corner (in CCW direction)
// -----------------------------------------------------------------------------
CPair CCorner::calculateF(double k, double c, double dx, CPair N, CPair relVel)
{
	double	dv = relVel*N;
	/* Kelvin-Voight */ return	-(k*dx-c*dv)*N;
	/* Kuwabara-Kono */ //return -(k*pow(dx,1.5) - c*sqrt(dx)*dv)*N;
	//return -(k*dx*exp(dx) - c*dv)*N;
}

// update the position and velocity
// --------------------------------
void CCorner::move(double dt)
{
	acc		= acc + CPair(0,g); 
	vel		= vel + acc*dt;
	pos		= pos + vel*dt;
	acc		= ZEROPAIR;	
	load	= 0.0;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Corner.h												//
//////////////////////////////////////////////////////////////////