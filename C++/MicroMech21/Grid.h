//////////////////////////////////////////////////////////////////
#ifndef GRID_H													//
#define GRID_H													//
// begin of Grid.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include <fstream.h>	// for file handling
#include <stdio.h>
#include <string>

extern double maxRadius;	// from main.cpp

extern double ceilY;		// from main.cpp
extern double inletFloorY;	// from main.cpp
extern double outletFloorY;	// from main.cpp
extern double inletX;		// from main.cpp
extern double tresholdX;	// from main.cpp
extern double outletX;		// from main.cpp
extern double vMean;		// from main.cpp
extern double minRadius;	// from main.cpp		


class CGrid
/*
 *	the Grid to perform Grid Search
 *  - space is divided into cells of sqrt(2)*Rmin width
 *	- normally one has to search within a range of 2*Rmax
 *  - since we have uniform particle size distribution Rmin=Rmax searchrange = 2 cells
 */
{	
public:	
	CGrid();
	~CGrid(){}
	void		update				( CBall ballA[], int nrBalls	);
	void		clearVelField		(								);
	void		addBall				( CBall ballA[], int ballNr		);
	void		draw				(								);
	CPairInt	cellNumber			( CPair pos						);
	CBall*		getBallP			( int gridCoX, int gridCoY		);
	int			getSearchRange		( double thisRadius, double maxRadius );
	void		drawCell			( int intX, int intY			);
	void		checkBoundaries		( CBall ballA[], int nrBalls	); 
	bool		isClear				( double xMin, double xMax		);
	void		imposeVelX			( double xMin, double xMax		);
	void		saveVelField		( const char* fileName			);
	int			getNrCells			(	){return nrCells.x*nrCells.y;}
private:
	

	double		cellWidth;
	CPairInt	nrCells;
	CBall***	cellA;			// 2D grid, filled with pointers to balls, NULL if empty
	CPair		lowerLeft;		// origin of the grid
	CPair		dim;			// the dimensions in x and y direction
	CPair**		sumVelA;
	int**		nrSamplesA;
};

// constructor
// -----------
CGrid::CGrid()
{
	// theoretically sqrt(2), 1.4 to be sure it's not too big
	cellWidth	= 1.4*minRadius; 
	double		margin = 5*getSearchRange(maxRadius,maxRadius)*cellWidth; // thus we make an extra border of 2 cells
	lowerLeft	= CPair(inletX-margin,outletFloorY-margin);
	dim			= CPair(outletX+margin,ceilY+margin)-lowerLeft;
	
	nrCells.x	= (int)ceil((dim.x)/cellWidth);
	nrCells.y	= (int)ceil((dim.y)/cellWidth);

	cellA = new CBall**[nrCells.x];
	sumVelA = new CPair*[nrCells.x];
	nrSamplesA = new int*[nrCells.x];
	for(int i=0; i<nrCells.x; ++i)
	{
		cellA[i] = new CBall*[nrCells.y];
		sumVelA[i] = new CPair[nrCells.y];
		nrSamplesA[i] = new int[nrCells.y];
		for(int j=0; j<nrCells.y; ++j)
		{
			cellA[i][j]	= NULL;
			sumVelA[i][j] = ZEROPAIR;
			nrSamplesA[i][j] = 0;
		}
	}
}

// calculate in which cell this point belongs
// ------------------------------------------
CPairInt CGrid::cellNumber(CPair pos)
{
	CPair nrPos = (pos-lowerLeft)/cellWidth;
	CPairInt cellNr;
	cellNr.x = (int)(int)floor(nrPos.x);
	cellNr.y = (int)(int)floor(nrPos.y);
	return cellNr;
}

// return the content of the grid at that position
// -----------------------------------------------
CBall* CGrid::getBallP(int gridCoX, int gridCoY)
{
	return cellA[gridCoX][gridCoY];
}

// return the range for contact detection
// --------------------------------------
int	CGrid::getSearchRange( double thisRadius, double maxRadius)
{
	return (int)ceil(thisRadius/cellWidth)+(int)ceil(maxRadius/cellWidth);
}

// add a ball to the grid
// ----------------------
void CGrid::addBall(CBall ballA[], int ballNr)
{
	CPairInt cellNr = CGrid::cellNumber(ballA[ballNr].getPos());
	cellA[cellNr.x][cellNr.y] = &ballA[ballNr];
}

// update the grid
// ---------------
void CGrid::update(CBall ballA[], int nrBalls)
{
	// A // clear everything: empty = NULL-pointer
	/////// --------------------------------------
	for(int i=0; i<nrCells.x; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			cellA[i][j]	= NULL;
		}
	}
	// B // insert the ball-pointers in the correct cell
	/////// --------------------------------------------
	for(i=0; i<nrBalls; ++i)
	{
		CPairInt cellNr = cellNumber(ballA[i].getPos());
		cellA[cellNr.x][cellNr.y] = &ballA[i];
		sumVelA[cellNr.x][cellNr.y] = sumVelA[cellNr.x][cellNr.y] + ballA[i].getVel();
		nrSamplesA[cellNr.x][cellNr.y] += 1;
	}
}

// clear the velocity field
// ------------------------
void CGrid::clearVelField()
{
	for(int i=0; i<nrCells.x; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			sumVelA[i][j]	= ZEROPAIR;
			nrSamplesA[i][j] = 0;
		}
	}
}

// save the velocity field
// -----------------------
void CGrid::saveVelField(const char* fileName)
{
	std::string stringNameX(fileName);
	std::string stringNameY(fileName);
	stringNameX.append("X.txt");
	stringNameY.append("Y.txt");
	//  char temp[33];
	//	itoa(i,temp,10);
	//	fileNameX.append(temp);
	const char* fileNameX = stringNameX.c_str();
	const char* fileNameY = stringNameY.c_str();

	ofstream outFileX;
	ofstream outFileY;
	outFileX.open(fileNameX, ios::out);
	outFileY.open(fileNameY, ios::out);
	if ((! outFileX)||(! outFileY)){cout << "Can't open output file.\n"; }
	
	CPair value;
	//for(int j=nrCells.y-1; j>=0; --j)
	for(int j=0; j<nrCells.y; ++j)
	{
		for(int i=0; i<nrCells.x; ++i)
		{
			if(nrSamplesA[i][j]>0) {
				value = sumVelA[i][j]/nrSamplesA[i][j];
				outFileX<<value.x<<" ";
				outFileY<<value.y<<" ";
			} else {
				outFileX<<0<<" ";
				outFileY<<0<<" ";
			}

		}
		outFileX<<endl;
		outFileY<<endl;
	}

/*	outFile	<<"# Average velocity in y-direction"<<endl
			<<"# --------------------------------"<<endl;
	for(i=0; i<nrCells.x; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			if(nrSamplesA[i][j]>0) {
				outFile<<sumVelA[i][j].y/nrSamplesA[i][j]<<" ";
			} else {
				outFile<<0<<" ";
			}
		}
		outFile<<endl;
	}
*/
}



// draw one cell
// -------------
void CGrid::drawCell(int intX, int intY)
{
	double posX = intX*cellWidth+lowerLeft.x;
	double posY = intY*cellWidth+lowerLeft.y;
	DrawRectangle(posX, posY, posX+cellWidth, posY+cellWidth); 
}

// draw the grid
// -------------
void CGrid::draw()
{
	// A // draw the contour of the entire grid
	/////// -----------------------------------
	CPair	  upperRight =  lowerLeft +   dim;
	DrawLine( lowerLeft.x,	lowerLeft.y,  upperRight.x, lowerLeft.y  );
	DrawLine( upperRight.x,	lowerLeft.y,  upperRight.x, upperRight.y );
	DrawLine( upperRight.x,	upperRight.y, lowerLeft.x,  upperRight.y );
	DrawLine( lowerLeft.x,	upperRight.y, lowerLeft.x,  lowerLeft.y  );

	// B // draw filled cells
	/////// -----------------
/*	for(int i=0; i<nrCells.x; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			if(cellA[i][j]==NULL) {
				// sit still and do nothing
			}
			else {
				//glColor3f(0,1,0);
				//drawCell(i,j); 
				glColor3f(0,0,1);
				cellA[i][j]->draw();
			}
		}
	}
	*/
}

// check the boundaries
// --------------------
void CGrid::checkBoundaries(CBall ballA[], int nrBalls)
{
	int boundaryThickness = 4;
	CPairInt cellNr = cellNumber(CPair(outletX,ceilY));
	for(int i=cellNr.x-boundaryThickness; i<cellNr.x+1; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			glColor3f(0,1,0);
			drawCell(i,j); 
			if(cellA[i][j]==NULL) {
				/* sit still and do nothing*/ }
			else {
				ballA[cellA[i][j]->nr].nr = -1;
				cellA[i][j] = NULL;
			}
		}
	}
}

// check the area between xMin and xMax is free of balls
// -----------------------------------------------------
bool CGrid::isClear(double xMin, double xMax)
{
	bool clear = true;
	//double	invRange = 1/(ceilY-inletFloorY);
	int cellNrMin = cellNumber(CPair(xMin,0)).x-1;
	int cellNrMax = cellNumber(CPair(xMax,0)).x+1;
	for(int i=cellNrMin; i<cellNrMax; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			if(cellA[i][j]!=NULL) {
				
				double posX = cellA[i][j]->getPos().x;
				

				if((posX >= xMin)&&(posX < xMax)) {
				//	double posY = cellA[i][j]->getPos().y;
				//	double s = (posY - inletFloorY)*invRange;
					double velX = cellA[i][j]->getVel().x;//1.2*vMean;//6.0*vMean*s*(1.0-s);
				cellA[i][j]->addF(CPair(GRAVITY,0));
					velX = 2*vMean;//1.01*fabs(velX);
					cellA[i][j]->setVelX(velX);
					cellA[i][j]->setVelY(0);
					glColor3f(1,0,0);
				drawCell(i,j);
				}
				if((posX >= xMin)&&(posX < xMax)) clear = false;
			}
		}
	}
	return clear;
}

// velocity filter: slow them down
// -------------------------------
void CGrid::imposeVelX(double xMin, double xMax)
{
	double	invRange = 1/(ceilY-inletFloorY);
	int cellNrMin = cellNumber(CPair(xMin,0)).x;
	int cellNrMax = cellNumber(CPair(xMax,0)).x+1;
	for(int i=cellNrMin; i<cellNrMax; ++i)
	{
		for(int j=0; j<nrCells.y; ++j)
		{
			if(cellA[i][j]!=NULL) {
				double posX = cellA[i][j]->getPos().x;
				if((posX > xMin)&&(posX < xMax)) {
					CPair pos = cellA[i][j]->getPos();
					double s = (pos.y - inletFloorY)*invRange;
					//double velX = 6.0*vMean*s*(1.0-s);
					double velX = cellA[i][j]->getVel().x;
					if(velX<6.0*vMean*s*(1.0-s)) { velX *= 1.001;}
					cellA[i][j]->setVelX(velX);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Grid.h												//
//////////////////////////////////////////////////////////////////