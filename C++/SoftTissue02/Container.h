//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "Corner.h"
#include "Cell.h"
#include "InputFile.h"

#include <fstream.h>	// for file handling

extern double floorHeight;	// from main.cpp
extern bool	  showForces;	// from main.cpp

class CContainer
/*
 * Contains an array of the corners which are linked to one another
 * each corner has its own dynamics which is determined by spring forces
 * connecting it to its two neighbours and a pressure force
 */
{	
public:	
	CContainer					( double newCornerMass, double newPressure	);
	~CContainer					(											){};
	void	draw				(											);
	void	move				( double dt									);
	void	calculateForces		( double k, double c, double pullForceX, double pullForceY);
	void	freeze				(											);
	void	pullUp();
//	void	save		();
private:
	CCell	*cellA;
	CCorner *cornerA;
	int		nrCells;
	int		nrCorners;
	double	ratioPxV;
	double	pressure;
	int		nrNorthFace;
	int		nrSouthFace;
	int		nrEastFace;
	int		nrWestFace;
	CCorner **northFace;
	CCorner **southFace;
	CCorner **eastFace;
	CCorner **westFace;
};

// to freeze the current configuration
// -----------------------------------
void CContainer::freeze()
{
	for(int i=0; i<nrCells; ++i)
	{
		cellA[i].freeze();
	}
}

// set the beginconditions
// -----------------------
CContainer::CContainer( double newCornerMass, double newPressure)
{
	int i;
	const char* nameIn = "data/inputFile.txt";
	InputFile iF(nameIn);
	// read the coordinates of the cornerpoints
	iF.setScope("-cornerpoints-");
	nrCorners = iF.getInt("nrcorners");
	cornerA = new CCorner[nrCorners];
	CPair	P;
	for(i=0; i<nrCorners; ++i)
	{
		P.x		= iF.getDouble();
		P.y		= iF.getDouble();
		cornerA[i].initialize(P,newCornerMass);
	}
	// read the configuration
	iF.setScope("-configuration-");
	nrCells = iF.getInt("nrcells");
	int *nrCornersA = new int[nrCells];
	int **cornerNrsA = new int*[nrCells];
	for(i=0; i<nrCells; ++i)
	{
		int nrCorners = iF.getInt("nrcorners");
		nrCornersA[i] = nrCorners;
		cornerNrsA[i] = new int[nrCorners];
		cornerNrsA[i][0] = iF.getInt("numbers");
		for(int j=1; j<nrCorners; ++j)
		{
			cornerNrsA[i][j] = iF.getInt();
		}
		iF.setScopeBegin(iF.location);
	}
	
	// read the northFace
	iF.setScope("-northface-");
	nrNorthFace = iF.getInt("nrcorners");
	northFace = new CCorner*[nrCorners];
	northFace[0] = &cornerA[iF.getInt("numbers")-1];
	for(i=1; i<nrNorthFace; ++i)
	{
		northFace[i] = &cornerA[iF.getInt()-1];
	}

	// read the southFace
	iF.setScope("-southface-");
	nrSouthFace = iF.getInt("nrcorners");
	southFace = new CCorner*[nrCorners];
	southFace[0] = &cornerA[iF.getInt("numbers")-1];
	for(i=1; i<nrSouthFace; ++i)
	{
		southFace[i] = &cornerA[iF.getInt()-1];
	}

	// read the eastFace
	iF.setScope("-eastface-");
	nrEastFace = iF.getInt("nrcorners");
	eastFace = new CCorner*[nrCorners];
	eastFace[0] = &cornerA[iF.getInt("numbers")-1];
	for(i=1; i<nrEastFace; ++i)
	{
		eastFace[i] = &cornerA[iF.getInt()-1];
	}

	// read the westFace
	iF.setScope("-westface-");
	nrWestFace = iF.getInt("nrcorners");
	westFace = new CCorner*[nrCorners];
	westFace[0] = &cornerA[iF.getInt("numbers")-1];
	for(i=1; i<nrWestFace; ++i)
	{
		westFace[i] = &cornerA[iF.getInt()-1];
	}



	// initialize all the cells with their respective corners
	cellA = new CCell[nrCells];
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].initialize(nrCornersA[i], cornerNrsA[i], cornerA, newPressure);
	}
	delete[] nrCornersA;
	delete[] cornerNrsA;
	
	// make this configuration the configuration in rest
	this->freeze();
}

/*void CContainer::save()
{
	const char* nameOut = "data/savedCornerPoints.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	if (! outFile){}
	for(int i=0; i<nrCorners; ++i)
	{

		outFile	<<cornerA[i].getPos().x<<" "
				<<cornerA[i].getPos().y<<endl;
	}
	outFile.close();
}


*/
// draw the current configuration
// ------------------------------
void CContainer::draw()
{	//draw the floor
//	glColor3f(0,0,1);
//	drawLine(CPair(0,floorHeight),CPair(1,floorHeight));

	for(int i=0; i<nrCells; ++i)
	{
		cellA[i].draw();
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
	// integrate the accelerations of the corners
	for(int i=0; i<nrCorners; ++i)
	{
		cornerA[i].move(dt);
	}
}

// calculate all the forces acting on the corners
// ----------------------------------------------
void CContainer::calculateForces(double k, double c, double pullForceX, double pullForceY)
/*
 * Calculate the force and modify the accelerations of the corners
 * N points from the begincorner to the endcorner (is CCW) and the
 * resulting F acts on the respective corners
 */
{
	// A // the corners are connected by springs
	/////// ------------------------------------
	int i;
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].calculateF(k,c);
	}

	CPair F;
	double dx; 
	CPair N;
	// B // contact with a floor on y = floorHeight
	/////// ---------------------------------------
	for(i=0; i<nrCorners; ++i)
	{
		dx = cornerA[i].getPos().y - floorHeight;
		if(dx<0)
		{
			N.x = 0; N.y = -1;
			F = cornerA[i].calculateF(k, c, dx, N, cornerA[i].getVel());
			cornerA[i].addF(-F);
		}
	}

	CPair midPoint;
	for(i=0; i<nrNorthFace; ++i)
	{
		glColor3f(1,1,0);
		F = pullForceY/nrNorthFace*CPair(0,1);
		if(showForces){ midPoint = northFace[i]->getPos();
		drawLine(midPoint,midPoint+0.005*F);}
		northFace[i]->addF(F);
	}
	for(i=0; i<nrEastFace; ++i)
	{
		glColor3f(1,1,0);
		F = pullForceX/nrEastFace*CPair(1,0);
		if(showForces){ midPoint = eastFace[i]->getPos();
		drawLine(midPoint,midPoint+0.005*F);}
		eastFace[i]->addF(F);
	}
	for(i=0; i<nrSouthFace; ++i)
	{
		glColor3f(1,1,0);
		F = pullForceY/nrSouthFace*CPair(0,-1);
		if(showForces){ midPoint = southFace[i]->getPos();
		drawLine(midPoint,midPoint+0.005*F);}
		southFace[i]->addF(F);
	}
	for(i=0; i<nrWestFace; ++i)
	{
		glColor3f(1,1,0);
		F = pullForceX/nrWestFace*CPair(-1,0);
		if(showForces){ midPoint = westFace[i]->getPos();
		drawLine(midPoint,midPoint+0.005*F);}
		westFace[i]->addF(F);
	}

	
	
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////