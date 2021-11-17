//////////////////////////////////////////////////////////////////
#ifndef CELL_H													//
#define CELL_H													//
// begin of Cell.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
//#include "InputFile.h"
#include "Spring.h"


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
	void	initialize(int newNr, int newNrCorners, int* cornerNrs, Corner* cornerA, 
							      int newNrSprings, int* springNrs, Spring* springA);
	void	draw();
	void	drawFull();	
	double	calculateSurface();
	void	save( ofstream outFile );
	void	save4SoftTissue( ofstream outFile, int nrActiveCorners );
	int		getNr(){ return nr; }
	int		getNrCorners(){return nrCorners;}
	int		getNrSprings(){return nrSprings;}
private:
	int		nr;
	Corner** pCornerA;
	Spring** pSpringA;
	int		nrCorners;
	int		nrSprings;
};

// set the initial values
// ----------------------
void Cell::initialize(int newNr, int newNrCorners, int* cornerNrs, Corner* cornerA, 
								 int newNrSprings, int* springNrs, Spring* springA)
{
	nr			= newNr;
	nrCorners	= newNrCorners;	
	nrSprings	= newNrSprings;
	pCornerA	= new Corner*[nrCorners];
	pSpringA	= new Spring*[nrSprings];

	int i;
	for(i=0; i<nrCorners; ++i){ 
		pCornerA[i] = &cornerA[cornerNrs[i]-1]; }
	for(i=0; i<nrSprings; ++i){
		pSpringA[i] = &springA[springNrs[i]-1]; }	
}

// draw the cell
// -------------
void Cell::draw() 
{ 
	int i;
	for( i=0; i<nrSprings; ++i ){ pSpringA[i]->draw(); }
}
void Cell::drawFull() 
{ 
	glBegin(GL_POLYGON);
		int i;
		for(i=0; i<nrCorners; ++i){
			glVertex2f(pCornerA[i]->getPos().x,pCornerA[i]->getPos().y); }
	glEnd();
}

// calculate the surface
// ---------------------
double Cell::calculateSurface()
{
	double	surface	= 0;	// the total surface
	double  IvectorProdI;
	Pair	A, B, C;
	A	  = pCornerA[0]->getPos();
	int i;
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

// save this cell configuration
// ----------------------------
void Cell::save4SoftTissue(ofstream outFile, int nrActiveCorners)
{
	outFile<<nr<<" "<<" nrcorners "<<2*nrCorners<<" numbers ";
	int i;
	for(i=0; i<nrCorners; ++i)
	{
		outFile<<pCornerA[i]->getNr4SoftTissue()<<" ";
	}
	for(i=0; i<nrCorners; ++i)
	{
		outFile<<pCornerA[i]->getNr4SoftTissue()+nrActiveCorners<<" ";
	}
	outFile<<endl;

}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Cell.h												//
//////////////////////////////////////////////////////////////////