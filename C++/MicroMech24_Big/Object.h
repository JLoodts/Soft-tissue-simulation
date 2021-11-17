//////////////////////////////////////////////////////////////////
#ifndef OBJECT_H												//
#define OBJECT_H												//
// begin of Object.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"

extern double time; // from main.cpp

class CObject
/*
 *	models an object constructed out of walls
 */
{	
public:	
	CObject(){}
	~CObject(){}
	void	initialize	( int newNr, double newRotFreq, CPair newPos	);
	void	draw		(								); 	
	void	calculateF	( double k,  double c, const int nrBalls, CBall *ballA);

	void	move		( /*double time	*/					);
	int		nr;				// nr in the objectArray of Container
private:
	CPair	pos;
	CPair	N;
	CPair	*wallPos;
	bool	*rotToRight;
	CWall	*wallA;
	int		nrWalls;
	CGrid	grid;
	int		maxNrWalls;
	double	vibrAmp;
	double	vibrFreq;
	double	rotFreq;
	double wallLength;
};

// set the initial values
// ----------------------
void CObject::initialize(int newNr, double newRotFreq, CPair newPos)
{
	nr			= newNr;
	nrWalls		= 0;
	maxNrWalls	= 2;
	pos			= newPos;
	wallLength = 0.01;

	vibrAmp		= 0.01;
	vibrFreq	= 30;
	rotFreq		= newRotFreq;

	wallA			= new CWall[maxNrWalls];
	wallPos	= new CPair[maxNrWalls];
	rotToRight		= new bool[maxNrWalls];

	CPair wallPosB, wallPosE;
	wallPos[nrWalls] = CPair(0.0,0.0);
	N = CPair(cos(rotFreq*time*2*PI+0.5*PI),sin(rotFreq*time*2*PI+0.5*PI));
	wallPosB = CPair(-N.x*wallLength,-N.y*wallLength) + pos + wallPos[nrWalls];
	wallPosE = CPair(N.x*wallLength,N.y*wallLength) + pos + wallPos[nrWalls];
	wallA[nrWalls].initialize(nrWalls, wallPosB,wallPosE);
	rotToRight[nrWalls] = true;
	nrWalls++;

	wallPos[nrWalls] = CPair(0.0,0.0);
	N = CPair(cos(rotFreq*time*2*PI),sin(rotFreq*time*2*PI));
	wallPosB = CPair(-N.x*wallLength,-N.y*wallLength) + pos + wallPos[nrWalls];
	wallPosE = CPair(N.x*wallLength,N.y*wallLength) + pos + wallPos[nrWalls];
	wallA[nrWalls].initialize(nrWalls, wallPosB,wallPosE);
	rotToRight[nrWalls] = true;
	nrWalls++;

	
}

void CObject::draw() 
{ 
	for(int i=0; i<nrWalls; ++i)
	{
		wallA[i].draw();
	}
}

// calculate the force 
// -------------------
void CObject::calculateF( double k,  double c, const int nrBalls, CBall *ballA)
{
	for(int i=0; i<nrWalls; ++i)
	{
		wallA[i].calculateF(k, c, nrBalls, ballA);
	}
}

// update the position and velocity
// --------------------------------
void CObject::move(/*double time*/)
{
	//CPair newPos = pos;
	CPair wallPosB, wallPosE;
	//newPos.x += vibrAmp*cos(vibrFreq*time*2*PI);

		N = CPair(cos(-rotFreq*time*2*PI+0.5*PI),sin(-rotFreq*time*2*PI+0.5*PI));
		wallPosB = CPair(-N.x*wallLength,-N.y*wallLength) + pos + wallPos[0];
		wallPosE = CPair(N.x*wallLength,N.y*wallLength) + pos + wallPos[0];
		wallA[0].move(wallPosB,wallPosE);

		N = CPair(cos(-rotFreq*time*2*PI),sin(-rotFreq*time*2*PI));
		wallPosB = CPair(-N.x*wallLength,-N.y*wallLength) + pos + wallPos[1];
		wallPosE = CPair(N.x*wallLength,N.y*wallLength) + pos + wallPos[1];
		wallA[1].move(wallPosB,wallPosE);
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Object.h												//
//////////////////////////////////////////////////////////////////