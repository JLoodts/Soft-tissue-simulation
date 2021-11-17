//////////////////////////////////////////////////////////////////
#ifndef CORNER_H												//
#define CORNER_H												//
// begin of Corner.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include <fstream.h>	// for file handling

extern double g;

/**
 *  The Corner class. 
 *	
 *	The most primitive element, all other elements like Springs and Cells contain Corners.
 *	Corners are the only things that have their own mass and dynamics
 */
class Corner
{	
public:
	Corner(){}
	~Corner(){}
	void	initialize( int newNr, CTriple newPos, double newMass );
	CTriple	getPos(){ return pos;	}
	CTriple	getVel(){ return vel; }
	CTriple	getF(){	return (acc + accP)*mass; }
	double	getMass(){ return mass; }
	void	setPosX( double newPosX ){ pos.x = newPosX; }
	void	setPosY( double newPosY ){ pos.y = newPosY; }
	void	setVel( CTriple newVel ){ vel = newVel; }
	void	setAcc( CTriple newAcc ){ acc = newAcc; }
	void	draw(){ drawPoint(pos); }
	void	move( double dt );
	void	addF( CTriple &F ){ acc = acc + (F/mass);	}
	void	addFtoAccP( CTriple &F){ accP = accP + (F/mass); }
	int		getNr(){ return nr; }
	void	save( ofstream outFile ){ outFile<<nr<<" "<<pos.x<<" "<<pos.y<<" "<<pos.z<<endl; }
	void	setBoundary(bool value){ boundary = value; }
private:
	int		nr;				/**< int Number of the corner. */
	CTriple	pos;			/**< CTriple Position [m] of the corner. */
	CTriple	vel;			/**< CTriple Velocity [m/s] of the corner. */
	CTriple	acc;			/**< CTriple Acceleration [m/s²] of the corner due to the spring force.*/
	CTriple	accP;			/**< CTriple Acceleration [m/s²] of the corner due to the internal pseudo pressure. */
	double	mass;			/**< double Mass [kg] of the corner. */
	bool	boundary;		// indicates that this point belongs to the boundary on which external forces are applied
};

// set the initial values
// ----------------------
void Corner::initialize(int newNr, CTriple newPos, double newMass)
{
	nr		= newNr;
	pos		= newPos;
	vel		= ZEROTRIPLE;
	acc		= ZEROTRIPLE;
	accP	= ZEROTRIPLE;
	mass	= newMass;
	boundary= false;
}


// update the position and velocity
// --------------------------------
void Corner::move(double dt)
{
	if(!boundary){
		acc		= acc + accP + CTriple(0,g,0); 
	} else {
		acc		= acc + CTriple(0,g,0);
	}
	vel		= vel + acc*dt;
	pos		= pos + vel*dt;
	acc		= ZEROTRIPLE;	
	accP	= ZEROTRIPLE;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Corner.h												//
//////////////////////////////////////////////////////////////////