//////////////////////////////////////////////////////////////////
#ifndef CELL_H													//
#define CELL_H													//
// begin of Cell.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "InputFile.h"
#include "Spring.h"

extern bool	showForces;		// from main.cpp
extern double scaleFactorForForces;	// from main.cpp
extern double dt;			// from main.cpp


#include <fstream.h>	// for file handling

/**
 * Cell class models a cell which is comprised of:
 * - an array of pointers to the cornerpoints
 * - an array of its sides, which are modeled as springs
 * - a measure for the internal pressure based on its surface/volume
 */
class Cell
{
public:	
	Cell(){}
	~Cell(){}
	void initialize(int newNr, int newNrCorners, int* cornerNrs, Corner* cornerA, 
							   int newNrSprings, int* springNrs, Spring* springA, const char* nameParam);
	void	drawFull();	
	void	calculateForce();
	double	calculateSurface();
	double	getSurface(){ return surfaceT; }
	double	getPressure(){ return pressure; }
	void	save( ofstream outFile );
	int		getNr(){ return nr; }

private:
	
	int		nr;
	double	k;				// the proportional constant to regulate a constant surface
	double	c;				// the damping constant to regulate a constant surface
	Corner** pCornerA;
	Spring** pSpringA;
	int		nrCorners;
	int		nrSprings;
	int		nrSpring2s;
	double	restSurface;		// the surface at the beginning wich is the equilibrium (target)
	double	surfaceTminusOne;	// the surface at the previous timestep
	double	surfaceT;			// the surface at the current timestep
	double	thickness;			// thickness of the tissue / cell (120 µm for onion)
	double	pressure;			// the current pressure in the cell (N/m²)
};

// set the initial values
// ----------------------
void Cell::initialize(int newNr, int newNrCorners, int* cornerNrs, Corner* cornerA, 
								 int newNrSprings, int* springNrs, Spring* springA, const char* nameParam)
{
	nr			= newNr;
	// put pointers to the corners in place
	nrCorners	= newNrCorners;
	pCornerA	= new Corner*[nrCorners];
	int i;
	for(i=0; i<nrCorners; ++i)
	{
		pCornerA[i] = &cornerA[cornerNrs[i]-1]; // -1 since the inputfile starts from 1 and c++ starts from 0
	}
	surfaceT = surfaceTminusOne = restSurface = this->calculateSurface();
//surfaceT = surfaceTminusOne = this->calculateSurface();
//restSurface = 1.01*surfaceT;	

	thickness = 120e-6; // for unicellular onion tissue
	pressure = 0;

	// put pointers to the springs in place
	nrSprings	= newNrSprings;
	pSpringA	= new Spring*[nrSprings];
	for(i=0; i<nrSprings; ++i)
	{
		pSpringA[i] = &springA[springNrs[i]-1]; // -1 since the inputfile starts from 1 and c++ starts from 0
	}

	{ // read from parameterFile.txt
//auto//		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameParam);

		iF.setScope("-parameters-");
		k	= iF.getDouble("kpressure");
		c	= 2*sqrt(k*cornerA[0].getMass()); 
		cout<<"Dpressure = "<<c<<endl;
	}
	

}


// calculate the force due to the spring force
// -------------------------------------------
void Cell::calculateForce()
{
	int i;
	double distance;
	Pair N, relVel, F, midPoint;

	// forces exerted by the internal pressure
	// ---------------------------------------
	surfaceTminusOne = surfaceT;
	surfaceT = this->calculateSurface();
	// 'surface velocity'/dt to make it [m²/s]
	double magnitude = k*(surfaceT-restSurface)+c*(surfaceT-surfaceTminusOne)/dt;
	pressure = magnitude/thickness;

	for(i=0; i<nrCorners-1; ++i)
	{
		N = pCornerA[i+1]->getPos() - pCornerA[i]->getPos();
		distance = VectorLength(N);
		N = N/distance; // unit length now
		F = Pair(-N.y,N.x); // rotation 90 deg to the left->pointing outward
		F = -(magnitude)*distance*F;

		if(showForces){// begin draw F
			Pair midPoint = pCornerA[i]->getPos()+N*0.5*distance;
			glColor3f(0,1,0);
			drawLine(midPoint ,midPoint+scaleFactorForForces*F);
		// end draw F
		}
		// devide the pressure force over the two corners
		pCornerA[i]->addFtoAccP(0.5*F);
		pCornerA[i+1]->addFtoAccP(0.5*F);
	}
	N = pCornerA[0]->getPos() - pCornerA[nrCorners-1]->getPos();
	distance = VectorLength(N);
	N = N/distance;
	F = Pair(-N.y,N.x); // rotation 90 deg to the left->pointing outward
	F = -(magnitude)*distance*F;

	if(showForces){// begin draw F
		midPoint = pCornerA[nrCorners-1]->getPos()+N*0.5*distance;
		glColor3f(0,1,0);
		drawLine(midPoint ,midPoint+scaleFactorForForces*F);
	// end draw F
	}
	// devide the pressure force over the two corners
	pCornerA[nrCorners-1]->addFtoAccP(0.5*F);
	pCornerA[0]->addFtoAccP(0.5*F);
		 
}


/*
void Cell::draw() 
{ 
	int i;
	for( i=0; i<nrSprings; ++i ){ pSpringA[i]->draw(); }
}*/

void Cell::drawFull() 
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
double Cell::calculateSurface()
{
	double	surface	= 0;	// the total surface
	int i;
	Pair	A, B, C;
	double IvectorProdI;
	A = pCornerA[0]->getPos();

	for(i=0; i<nrCorners-1; ++i)
	{
		B = pCornerA[i]->getPos()-A;
		C = pCornerA[i+1]->getPos()-A;
		IvectorProdI = B.x*C.y-B.y*C.x;

		surface += 0.5*IvectorProdI;
	}

	return -surface; // minus since we went clockwise
}

// save this cell configuration
// ----------------------------
void Cell::save(ofstream outFile)
{
	outFile<<nr<<" "<<" nrcorners "<<nrCorners<<" numbers ";
	int i;
	for(i=0; i<nrCorners; ++i)
	{
		outFile<<pCornerA[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"  "<<" nrsprings "<<nrSprings<<" numbers ";
	for(i=0; i<nrSprings; ++i)
	{
		outFile<<pSpringA[i]->getNr()<<" ";
	}
	outFile<<endl;
}




//////////////////////////////////////////////////////////////////
#endif															//
// end of Cell.h												//
//////////////////////////////////////////////////////////////////