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
	void	initialize( int newNr, Pair newPos, double newMass );
	Pair	getPos(){ return pos;	}
	Pair	getVel(){ return vel; }
	double	getMass(){ return mass; }
	void	setPosX( double newPosX ){ pos.x = newPosX; }
	void	setVel( Pair newVel ){ vel = newVel; }
	void	setAcc( Pair newAcc ){ acc = newAcc; }
	void	draw(){ drawPoint(pos); }
	void	move( double dt );
	void	addF( Pair &F ){ acc = acc + (F/mass);	}
	void	addFtoAccP( Pair &F){ accP = accP + (F/mass); }
	int		getNr(){ return nr; }
	void	save( ofstream outFile ){ outFile<<nr<<" "<<pos.x<<" "<<pos.y<<endl; }
	void	setBoundary(bool value){ boundary = value; }
private:
	int		nr;				/**< int Number of the corner. */
	Pair	pos;			/**< Pair Position [m] of the corner. */
	Pair	vel;			/**< Pair Velocity [m/s] of the corner. */
	Pair	acc;			/**< Pair Acceleration [m/s²] of the corner due to the spring force.*/
	Pair	accP;			/**< Pair Acceleration [m/s²] of the corner due to the internal pseudo pressure. */
	double	mass;			/**< double Mass [kg] of the corner. */
	bool	boundary;		// indicates that this point belongs to the boundary on which external forces are applied
};

// set the initial values
// ----------------------
void Corner::initialize(int newNr, Pair newPos, double newMass)
{
	nr		= newNr;
	pos		= newPos;
	vel		= ZEROPAIR;
	acc		= ZEROPAIR;
	accP	= ZEROPAIR;
	mass	= newMass;
	boundary= false;
}

// update the position and velocity
// --------------------------------
void Corner::move(double dt)
{
	if(!boundary){
		acc		= acc + accP + Pair(0,g); 
	} else {
		acc		= acc + Pair(0,g);
	}
	vel		= vel + acc*dt;
	pos		= pos + vel*dt;
	acc		= ZEROPAIR;	
	accP	= ZEROPAIR;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Corner.h												//
//////////////////////////////////////////////////////////////////