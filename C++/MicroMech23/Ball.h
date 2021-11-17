//////////////////////////////////////////////////////////////////
#ifndef BALL_H													//
#define BALL_H													//
// begin of Ball.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"

//extern double ballMass; // from Main.h
extern double swarmImpact; // from Main.h
extern CPair  gravity;	// from Main.h
extern double minRadius;	// minimum radius of the balls
extern double maxRadius; // radius of the balls
extern double minMass;
extern double maxMass;

class CBall
/*
 *	models a ball that constitutes the fluid:
 *  - soft collisions with a spring repulsive force
 *	  in response to the virtual overlap in the normal direction
 *	- these are the only things that move
 */
{	
public:	
	CBall(){}
	~CBall(){}
	void	initialize	( int newNr, CPair newPos, CPair newVel, double newRadius, double newMass	);
//	void	redirect	( CPair swarmDirection			);
	CPair	getPos		(								) { return pos;				}
	CPair	getVel		(								) { return vel;				}
	double	getRadius	(								) { return radius;			}
//	void	bounce		( CBall *pPartner );
	CPair	calculateF	( double k,  double c, double distance, CPair N, CPair relVel );
	void	setVelX		( double velX					 ) {vel.x = velX;}
	void	setVelY		( double velY					) {vel.y = velY;}
	void	draw		(								); 
	void	draw		(	double fraction				) { DrawDisc(pos,fraction*radius);	}
	void	drawPoint	(								) { DrawPoint(pos);			}
	CPair	move		( double dt						);
	void	addF		( CPair &F						) { acc = acc + (F/m);	}
	int		nr;				// nr in the ballArray of Container
private:
	double	radius;			// radius
	float	color;			// between 0 and 1
	CPair	pos;			// position
	CPair	vel;			// velocity
	CPair	acc;			// acceleration
	double	m;				// mass

};

// set the initial values
// ----------------------
void CBall::initialize(int newNr, CPair newPos, CPair newVel, double newRadius, double newMass)
{
	nr		= newNr;
	pos		= newPos;
	vel		= newVel;
	acc		= ZEROPAIR;
	radius	= newRadius;
	m		= newMass;
	color	= (radius - minRadius)/(maxRadius-minRadius);
	//color	= (m - minMass)/(maxMass-minMass);
}

extern double velMax;
void CBall::draw() 
{ 
	/* for the brazilian nut
	if(nr!=0) {	
		color	= (velMax + vel.y)/(2.0*velMax);
		glColor3f(color,0,1-color);
	} else {glColor3f(1,0.3,0);}
	*/
	color	= (velMax + vel.x)/(2.0*velMax);
	glColor3f(color,0,1-color);	
	DrawDisc(pos,radius);	
}

// redirect the velocity by the swarm direction
// --------------------------------------------
/*void CBall::redirect(CPair swarmDirection)
{
	double IIvelII = VectorLength(vel);
	if(IIvelII==0){vel = ZEROPAIR;}
	else{vel = vel/IIvelII;}
	vel.x = (1.0 - swarmImpact) * vel.x + swarmImpact * swarmDirection.x;
	vel.y = (1.0 - swarmImpact) * vel.y + swarmImpact * swarmDirection.y;
//	acc.x = swarmDirection.x;
//	acc.y = swarmDirection.y;
	double IIvelNewII = VectorLength(vel);
	if(IIvelNewII==0){vel = ZEROPAIR;}
	else{	vel = IIvelII*(vel/IIvelNewII);}
}

extern double c;
void CBall::bounce(CBall *pPartner)
{
	CPair N = pPartner->getPos() - pos;
	CPair velN = N*(N*vel);
	CPair velT = vel - velN;
	CPair velNnew = N*(N*pPartner->getVel());
	vel = (velNnew+velT);
}
*/
// calculate the force due to the walls of other balls
// ---------------------------------------------------
CPair CBall::calculateF(double k, double c, double overlap, CPair N, CPair relVel)
{
	double	dx = overlap;
	double	dv = relVel*N;
	return	-(k*dx-c*dv)*N;
}

// update the position and velocity
// --------------------------------
CPair CBall::move(double dt)
{
	vel		= vel + (acc+gravity)*dt;
//	cout<<"velocity"<<VectorLength(vel)<<endl;
	pos		= pos + vel*dt;
	acc		= ZEROPAIR;	
	return vel;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Ball.h												//
//////////////////////////////////////////////////////////////////