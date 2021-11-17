//////////////////////////////////////////////////////////////////
#ifndef WALL_H													//
#define WALL_H													//
// begin of Wall.h												//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "Grid.h"

#include <list>
typedef std::list<CPairInt> LIST;

class CWall
/*
 *	models a wall:
 *  - begin and endpoint determine the wall
 *	- are fixed
 */
{	
public:	
	CWall(){}
	~CWall(){}
	void	initialize	( int newNr, CPair beginPos, CPair endPos, CGrid* newGrid );
	void	initialize	(int newNr, CPair newBeginPos, CPair newEndPos);
	void	calculateF	( double k,  double c, const int nrBalls, CBall *ballA);
	void	draw();
	void	move(CPair newBeginPos, CPair newEndPos);
	int		nr;				// nr in the wallArray of Container
private:
	CPair	beginPos;		// begin position
	CPair	endPos;			// end position
	LIST	gridCell;		// list of all the gridCells in the neighbourhood of the wall
	LIST::iterator listIt;	// iterator to traverse the above list
	CPair	W;				// unit vector pointing from begin to end
	double	length;			// length of the wall
	CGrid*	grid;
	bool	fixed;			// if fixed the not all the grid must be checked
};

void CWall::initialize(int newNr, CPair newBeginPos, CPair newEndPos)
{
	nr = newNr;
	beginPos = newBeginPos;
	endPos = newEndPos;
	W = endPos - beginPos;
	length = VectorLength(W);
	W = W/length;
	fixed = false; //because no grid is passed
}

void CWall::initialize(int newNr, CPair newBeginPos, CPair newEndPos, CGrid* newGrid)
{
	this->initialize(newNr, newBeginPos, newEndPos);

	fixed = true; // otherwise the grid would not have been passed

	grid = newGrid;
	// fill the list with the appropriate grid cells
	double stepSize = grid->getCellWidth();
	int nrSteps = (int)(length/stepSize);
	CPair pointOnWall;
	int border = 2;
	CPairInt cellNumber;
	for(int i=0; i<nrSteps; ++i)
	{
		pointOnWall = beginPos + i*stepSize*W;
		cellNumber = grid->cellNumber(pointOnWall);
		for(int j=-border; j<border; ++j)
		{
			for(int k=-border; k<border; ++k)
			{
// TO DO //		// still need to remove duplicate entries!!!
				gridCell.push_back(CPairInt(cellNumber.x+j, cellNumber.y+k));
			}
		}
		
	}
	
}

void CWall::calculateF(double k,  double c, const int nrBalls, CBall *ballA)
{
	CPair newP, closestP;
	double l;
	if(!fixed) {
		// check every ball with this wall
		for(int i=0; i<nrBalls; ++i)
		{
			newP = ballA[i].getPos() - beginPos;
			l = newP.x*W.x + newP.y*W.y;
			if(l<=0){ closestP = beginPos;}
			else	{ if(l >= length) {	closestP = endPos;}
					  else		  {	closestP = l*W + beginPos;}
					}
			// calculate the force and modify the accelerations
			CPair N = ballA[i].getPos() - closestP;
			double distance = VectorLength(N);
			double overlap = distance - ballA[i].getRadius();
			if(overlap<0)
			{
				// for shaking box
					extern double amplitude;
					extern double frequency;
					extern double time;
					double wallVel = -amplitude*frequency*2*PI*sin(frequency*time*2*PI);
				N = N/distance;
				CPair F = ballA[i].calculateF(k, c, overlap, N, ZEROPAIR/*CPair(0,wallVel)-ballA[i].getVel()*/);
				ballA[i].addF(F);
			}
		}
	} else {
		// look only in the neighbourhood
		CBall *ball;
		for (listIt=gridCell.begin(); listIt != gridCell.end(); ++listIt)  
		{
			ball = grid->getBallP(listIt->x,listIt->y);
			if(ball==NULL) {
				// sit still and do nothing
			}
			else {
				newP = ball->getPos() - beginPos;
				l = newP.x*W.x + newP.y*W.y;
				if(l<=0){ closestP = beginPos;}
				else	{ if(l >= length) {	closestP = endPos;}
						  else		  {	closestP = l*W + beginPos;}
						}
				// calculate the force and modify the accelerations
				CPair N = ball->getPos() - closestP;
				double distance = VectorLength(N);
				double overlap = distance - ball->getRadius();
				if(overlap<0)
				{
					
					N = N/distance;
					CPair F = ball->calculateF(k, c, overlap, N, -ball->getVel());
					ball->addF(F);
				}
			}
		}
	}
}

void CWall::draw()
{
	DrawLine(beginPos, endPos);
/*	for (listIt=gridCell.begin(); listIt != gridCell.end(); ++listIt)  
	{
		glColor3f(1.0,0.0,0.0);
		grid->drawCell(listIt->x, listIt->y);
	}
*/
}

void CWall::move(CPair newBeginPos, CPair newEndPos)
{
	beginPos = newBeginPos; 
	endPos = newEndPos;

	W = endPos - beginPos;
	length = VectorLength(W);
	W = W/length;
	fixed = false;
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Wall.h												//
//////////////////////////////////////////////////////////////////