//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "Main.h"			// for Random()
#include "Ball.h"
#include "Wall.h"
#include "Object.h"
#include "Grid.h"

extern double ceilY;		// from main.cpp
extern double inletFloorY;	// from main.cpp
extern double outletFloorY;	// from main.cpp
extern double inletX;		// from main.cpp
extern double tresholdX;	// from main.cpp
extern double outletX;		// from main.cpp
extern double vMean;		// from main.cpp
extern double friction;		// from main.cpp
//extern double dt;			// from main.cpp
//extern double minMass;
//extern double maxMass;


class CContainer
/*
 * Contains an array of the balls and walls. It can calculate all 
 * interaction forces
 */
{	
public:	
	CContainer					(						){}
	CContainer					( int newMaxNrBalls, int newMaxNrWalls, int newMaxNrObjects, double newMinRadius, double newMaxRadius,double newMinMass, double newMaxMass  );
	~CContainer					(						){}
	void	initialize			(						);
	void	draw				(						);
	void	clearVelField		(){grid.clearVelField();};
	void	saveVelField		(const char* fileName){grid.saveVelField(fileName);};
	double	move				( double dt				);
	void	calculateForces		( double k, double c	);
	void	save				(						);
	double	calculatePressure	(						);
	double	calculateSurface	(						);
	void	deleteBall			( int ballNr			);
	void	createBall			( double x, double y	);
	void	checkBoundaries		(						);
	int		getNrBalls			(){return nrBalls;}
	int		getNrCells			(){return grid.getNrCells();}
private:
	CBall	*ballA;
	CWall	*wallA;
	CObject *objectA;
	CGrid	grid;
	int		nrBalls;
	int		nrWalls;
	int		nrObjects;
	int		maxNrBalls;
	int		maxNrWalls;
	int		maxNrObjects;
	double	minRadius;
	double	maxRadius;
	double	minMass;
	double	maxMass;
	void	fillBox(double xMin, double xMax, double yMin, double yMax);
};

// constructor, initialize() must be called as well!!!
// ---------------------------------------------------
CContainer::CContainer(int newMaxNrBalls, int newMaxNrWalls, int newMaxNrObjects, double newMinRadius, double newMaxRadius,double newMinMass, double newMaxMass)
{
	nrBalls		= 0;
	nrWalls		= 0;
	nrObjects	= 0;
	maxNrBalls	= newMaxNrBalls;
	maxNrWalls	= newMaxNrWalls;
	maxNrObjects = newMaxNrObjects;
	minRadius	= newMinRadius;
	maxRadius	= newMaxRadius;
	minMass	= newMinMass;
	maxMass	= newMaxMass;
	ballA		= new CBall[maxNrBalls];
	wallA		= new CWall[maxNrWalls];
	objectA		= new CObject[maxNrObjects];
}

void CContainer::fillBox(double xMin, double xMax, double yMin, double yMax)
{
	double	invRange = 1/(yMax-yMin);
	double	beginPosX = xMin + 1.001*maxRadius;
	double	beginPosY = yMin + 1.001*maxRadius; /*for the big nut */// + 8.001*maxRadius;
	double	beginVelX = 0.0;
	double	beginVelY = 0.0;
	double	radius	  = Random(minRadius, maxRadius);
	double	mass	  = Random(minMass, maxMass);
	bool	stop	  = false;
	
	// fill the inlet
	
	// one big Brazilian nut
	/*	double bigRadius = 4*maxRadius;
	CPair bigBeginPos = CPair(0.5*(xMin+xMax),yMin + 1.001*bigRadius);
	double bigMass = 5*mass;//*pow(bigRadius/radius,3.0);
	ballA[nrBalls].initialize(nrBalls, bigBeginPos, CPair(beginVelX,beginVelY),
		bigRadius,bigMass);
			nrBalls++;
	// the rest of the balls
	*/
	if(nrBalls >= maxNrBalls){stop = true;}
	while((beginPosY<yMax-maxRadius)&&(!stop))
	{
		while((beginPosX<xMax-maxRadius)&&(!stop))
		{
			//double s = (beginPosY - inletFloorY)*invRange;
			//beginVelX = 6.0*vMean*s*(1.0-s);
			beginVelX = Random(-vMean,vMean);
			beginVelY = Random(-vMean,vMean);
			ballA[nrBalls].initialize(nrBalls, CPair(beginPosX,beginPosY),
									  CPair(beginVelX,beginVelY),radius,mass);
			nrBalls++;
			if(nrBalls >= maxNrBalls){stop = true;}
			beginPosX += maxRadius + 1.001*maxRadius;
			radius = Random(minRadius, maxRadius);
			//mass	  = pow(radius/maxRadius,3.0)*ballMass;
			mass = Random(minMass,maxMass);
		}
		beginPosY += maxRadius + 1.001*maxRadius;	
		beginPosX = xMin + 1.001*maxRadius;
	}
}
// set the beginconditions
// -----------------------
void CContainer::initialize()
/*
 * Calculate beginpositions of the balls
 * at the inlet there is a parabolic velocity profile with vMean
 * as the mean velocity
 */
{
	// the inlet
	fillBox(0, 0.1, 0.05, 0.07);
	fillBox(0, 0.056, 0.0, 0.025);

	cout<<"Number of balls: "<<nrBalls<<endl;

	grid.update(ballA, nrBalls);

	// initialize the walls
	extern double rightPos;
	wallA[nrWalls].initialize(nrWalls,CPair(0.0,0.05),CPair(0.05-rightPos, 0.05), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.05,0.05),CPair(0.1-3.0*rightPos, 0.05), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.05,0.05),CPair(0.04, 0.04), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.1,0.05),CPair(0.03, 0.03), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.015,0.025),CPair(0.02, 0.03), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.02, 0.03),CPair(0.025, 0.025), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.025, 0.025),CPair(0.015, 0.025), &grid);
	nrWalls++;



	wallA[nrWalls].initialize(nrWalls,CPair(0,0.07),CPair(0.1,0.07), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.1,0.07),CPair(0.1,0), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0.1,0),CPair(0,0), &grid);
	nrWalls++;
	wallA[nrWalls].initialize(nrWalls,CPair(0,0),CPair(0,0.07), &grid);
	nrWalls++;
	



	objectA[nrObjects].initialize(nrObjects, -3, CPair(0.066, 0.011));
	nrObjects++;
	objectA[nrObjects].initialize(nrObjects, 3, CPair(0.084, 0.011));
	nrObjects++;
/**///	objectA[nrObjects].initialize(nrObjects, 1, CPair(0.0, 0.1));
//	nrObjects++;



}

// draw the current configuration
// ------------------------------
void CContainer::draw()
{	
	int i;
	// A // the grid in green
	/////// -----------------
	glColor3f(0,1,0);
	grid.draw();

	// B // the walls in blue
	/////// -----------------
	glColor3f(0,0,1);
	for(i=0; i<nrWalls; ++i)
	{
		wallA[i].draw();
	}

	// C // the balls in blue
	/////// -----------------
	glColor3f(0,0,1);
	for(i=0; i<nrBalls; ++i)
	{
		ballA[i].draw();
	}

	// D // the objects in red
	/////// -----------------
	glColor3f(1,0,0);
	for(i=0; i<nrObjects; ++i)
	{
		objectA[i].draw();
	}

// keep track of ball nr 50
//	glColor3f(0.5,0.5,0.5);
//	ballA[50].draw(ballRadius);

}

// move the balls by integrating the equation of motion
// ------------------------------------------------------
double CContainer::move(double dt)
/*
 * The balls updated their accelerations in the calculateForces step
 * now their positions will be updated
 */
{
	CPair vel(0,0);
	double velMax = 0;
	for(int i=0; i<nrBalls; ++i)
	{
		vel = ballA[i].move(dt);
		if(fabs(vel.x)>velMax){ 
			velMax = fabs(vel.x);}
		if(fabs(vel.y)>velMax){ 
			velMax = fabs(vel.y);}
	}
	
	grid.update(ballA,nrBalls);
	for(i=0; i<nrObjects; ++i)
	{
		objectA[i].move();
	}

	extern double rightPos;
	wallA[0].move(CPair(0.0,0.05),CPair(0.05-rightPos, 0.05));
	wallA[1].move(CPair(0.05,0.05),CPair(0.1-3.0*rightPos, 0.05));
	
//	checkBoundaries();
//	grid.imposeVelX(tresholdX-4.1*maxRadius,tresholdX);
//	grid.imposeVelX(inletX+5*maxRadius,tresholdX-2*maxRadius);
	return velMax;
}

// calculate all the forces acting on the corners
// ----------------------------------------------
void CContainer::calculateForces(double k, double c)
/*
 * Calculate the force and modify the accelerations of the balls
 */
{
	// A // the balls interact by soft-collisions
	/////// -------------------------------------
	CPair N, F, relVel;
	double distance, overlap;
	int i;
	// check every possible contact
/*	for(i=1; i<nrBalls; ++i)
	{
		for(int j=0; j<i; ++j)
		{
			if((i!=0)||(j!=0)) 
			{
				N = ballA[j].getPos() - ballA[i].getPos(); // points from i to j
				distance = VectorLength(N);
				overlap = -(distance - ballA[i].getRadius() - ballA[j].getRadius());
				if(overlap>=0) { 
	//				glColor3f(0,1,0); DrawPoint(ballA[i].getPos()+0.5*N);
					N = N/distance; // unit length now
					relVel = ballA[j].getVel() - ballA[i].getVel();
					F = ballA[i].calculateF(k, c, overlap, N, relVel);
					ballA[i].addF(F); ballA[j].addF(-F);
				}
			}
		}
	}
*/
	// grid search
	CBall* pPartner;
	for(int n=0; n<nrBalls; ++n)
	{
//		CPair swarmVel(0,0);
//		double IIswarmVelII;
		CPair partnerVel;
		int searchRange = 1*grid.getSearchRange(ballA[n].getRadius(),maxRadius);
		CPairInt gridCo = grid.cellNumber(ballA[n].getPos());
		for(int i=-searchRange; i<=+searchRange; ++i)
		{
			for(int j=-searchRange; j<=+searchRange; ++j)
			{
				
				pPartner = grid.getBallP(gridCo.x+i,gridCo.y+j);
				//if(n==1){glColor3f(0.8,0.8,0.8); grid.drawCell(gridCo.x+i,gridCo.y+j);}
				if(pPartner!=NULL) 
				{
					//if(n==1) {glColor3f(1,0,0);pPartner->draw(2);}
					if((i!=0)||(j!=0)) {
						partnerVel = pPartner->getVel();
//						swarmVel = swarmVel + partnerVel;
						N = pPartner->getPos() - ballA[n].getPos(); // points from n to j
						distance = VectorLength(N);
						double sumRadii = ballA[n].getRadius() + pPartner->getRadius();
						overlap = -(distance - sumRadii);
						if(overlap>0) { 
							N = N/distance; // unit length now
							relVel = partnerVel - ballA[n].getVel();
							F = ballA[n].calculateF(k, c, overlap, N, relVel);
							ballA[n].addF(F); pPartner->addF(-F);
						}/*else {
							if((distance>=sumRadii)&&(distance<sumRadii+1.1*maxRadius)) { 
//								glColor3f(0,1,0); DrawPoint(ballA[n].getPos()+0.5*N);
								N = N/distance; // unit length now
								relVel = partnerVel - ballA[n].getVel();
								F = -0.00005*ballA[n].calculateF(k, c, distance, N, ZEROPAIR);
								ballA[n].addF(F); pPartner->addF(-F);
							}
						}*/
					}
				}
			}
		}
	/*	IIswarmVelII = VectorLength(swarmVel);
		if(IIswarmVelII!=0){ swarmVel = swarmVel/IIswarmVelII;}
		else{swarmVel = CPair(0.70710678,0.70710678);}
		ballA[n].redirect(swarmVel);
	*/
	}
	
	// B // contact with a wall
	/////// -------------------
	for(i=0; i<nrWalls; ++i)
	{
		wallA[i].calculateF(k, c, nrBalls, ballA);
	}

/*	for(i=0; i<nrBalls; ++i)
	{
		
		double overlap;
		double x = ballA[i].getPos().x;
		double y = ballA[i].getPos().y;
		double radius = ballA[i].getRadius();

		// contact with upper wall
		// -----------------------
		if(y >= ceilY-radius) {
			overlap = y + radius - ceilY ;
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawPoint(x,ceilY);
				N.x = 0; N.y = 1;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
				F.x = 0;
				ballA[i].addF(F);
			//	if(x>inletX+3*maxRadius)
			//	{ballA[i].setVelX(friction*ballA[i].getVel().x);}
			//	ballA[i].bounce();
			}
		}

		// contact with lower left wall
		// ----------------------------
		if((x < tresholdX)&&(y >= inletFloorY)) {
			overlap = inletFloorY - y + radius;
			if(overlap>=0) { // contact with lower left wall
//				glColor3f(1,0,0); DrawPoint(x,inletFloorY);
				N.x = 0; N.y = -1;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
			//	F.x = 0;
				ballA[i].addF(F);
				{ballA[i].setVelX(friction*ballA[i].getVel().x);
				ballA[i].setVelY(friction*ballA[i].getVel().y);}
			}
		} 
		
		// contact with treshold wall
		// --------------------------
		if(y < inletFloorY) {
			overlap = tresholdX - x + radius;
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawPoint(tresholdX,ballA[i].getPos().y);
				N.x = -1; N.y = 0;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
			//	F.y = 0;	
				ballA[i].addF(F);
				{ballA[i].setVelY(friction*ballA[i].getVel().x);}
			}
		}

		// contact with corner point treshold wall
		// ---------------------------------------
		if(((y >= inletFloorY)&&(y <= inletFloorY + radius))&&
		   ((x >= tresholdX  )&&(x <= tresholdX   + radius))) {
			overlap = -(Distance(x,y,tresholdX,inletFloorY) - radius);
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawDisc(CPair(tresholdX,inletFloorY),0.01);
				N.x = -(x-tresholdX); N.y = -(y-inletFloorY);
				N = N/VectorLength(N);
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
				ballA[i].addF(F);
			}
		}

		// contact with lower right wall
		// -----------------------------
		if(y <= outletFloorY+radius) {		
			overlap = outletFloorY - y + radius;
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawPoint(ballA[i].getPos().x,outletFloorY);
				N.x = 0; N.y = -1;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
				F.x = 0;
				ballA[i].addF(F);
		
			//	if(x>inletX+3*maxRadius)
			//	{ballA[i].setVelX(friction*ballA[i].getVel().x);}
			//	ballA[i].bounce();
			}
		}
		// contact with inlet wall
		// -----------------------
		if(x <= inletX+radius) {
			overlap = inletX - x +  radius;
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawPoint(inletX,ballA[i].getPos().y);
				N.x = -1; N.y = 0;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
				F.y = 0;
				ballA[i].addF(F);
			//	{ballA[i].setVelX(friction*ballA[i].getVel().x);}
			//	ballA[i].bounce();
			}
		}

		// contact with outlet wall
		// ------------------------
		if(x >= outletX-radius) {
			overlap =  x +  radius - outletX;
			if(overlap>=0) { 
//				glColor3f(1,0,0); DrawPoint(outletX,ballA[i].getPos().y);
				N.x = 1; N.y = 0;
				F = ballA[i].calculateF(k, c, overlap, N, -ballA[i].getVel());
				F.y = 0;
				ballA[i].addF(F);
			//	{ballA[i].setVelX(friction*ballA[i].getVel().x);}
			//	ballA[i].bounce();
			}
		}

	
	}
*/
	// C // contact with an object
	/////// ----------------------
	for(i=0; i<nrObjects; ++i)
	{
		objectA[i].calculateF(k, c, nrBalls, ballA);
	}

}

// save the current configuration
// ------------------------------
void CContainer::save()
{}

// delete a ball and switch its position
// -------------------------------------
void CContainer::deleteBall( int ballNr	)
{
	CBall* temp;
	temp = &ballA[nrBalls-1];
	ballA[ballNr] = *temp;
	ballA[ballNr].nr = ballNr;
	nrBalls--;
}

// create a ball 
// -------------
void CContainer::createBall(double x, double y)
{
	if(nrBalls<maxNrBalls)
	{
		CPair pos(x,y);
		CPair vel(Random(0,vMean),0);
		double radius = Random(minRadius, maxRadius);
		double mass	  = Random(minMass, maxMass);
		ballA[nrBalls].initialize(nrBalls,pos,vel,radius,mass);
/* niet echt nodig? */	//	grid.addBall(ballA,nrBalls);
		
		nrBalls++;
	}
}

// delete or create balls
// -------------------------------------
void CContainer::checkBoundaries()
{
	grid.checkBoundaries(ballA, nrBalls);
	int i;
	for( i=0; i<nrBalls; ++i)
	{
		if(ballA[i].nr==-1)
		{
			deleteBall(i);
//			cout<<"deleted: "<<i<<" nrBalls: "<<nrBalls<<endl;
				
		}
	}

	if(grid.isClear(inletX-maxRadius,inletX+3.003*maxRadius)) {
		bool	stop	  = false;
		double beginPosY = inletFloorY + 1.001*maxRadius;
		double radius = Random(minRadius, maxRadius);
		while((beginPosY<ceilY-maxRadius)&&(!stop))
		{
			createBall(inletX + 1.001*radius,beginPosY);
			
//			cout<<"created: "<<nrBalls-1<<" nrBalls: "<<nrBalls<<endl;
			beginPosY += 1.001*radius + radius;
			if(nrBalls == maxNrBalls){stop = true;}
			radius = Random(minRadius, maxRadius);
		}
	}

}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////