//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "Main.h"			// for Random()
#include "Ball.h"
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
	CContainer					( int newMaxNrBalls, double newMinRadius, double newMaxRadius, double newMinMass, double newMaxMass  );
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
	CGrid	grid;
	int		nrBalls;
	int		maxNrBalls;
	double	minRadius;
	double	maxRadius;
	double	minMass;
	double	maxMass;
};

// constructor, initialize() must be called as well!!!
// ---------------------------------------------------
CContainer::CContainer(int newMaxNrBalls, double newMinRadius, double newMaxRadius,double newMinMass, double newMaxMass)
{
	nrBalls		= 0;
	maxNrBalls	= newMaxNrBalls;
	minRadius	= newMinRadius;
	maxRadius	= newMaxRadius;
	minMass	= newMinMass;
	maxMass	= newMaxMass;
	ballA		= new CBall[maxNrBalls];
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
	double	invRange = 1/(ceilY-inletFloorY);
	double	beginPosX = inletX + 1.001*maxRadius;
	//double	beginPosY = inletFloorY + 1.001*ballRadius;
	double	beginPosY = inletFloorY + 1.001*maxRadius;
	double	beginVelX = 0.0;
	double	beginVelY = 0.0;
	double	radius	  = Random(minRadius, maxRadius);
	double	mass	  = Random(minMass, maxMass);
	bool	stop	  = false;
	
	// fill the inlet
	while((beginPosY<ceilY-maxRadius)&&(!stop))
	{
		while((beginPosX<tresholdX-maxRadius)&&(!stop))
		{
			//double s = (beginPosY - inletFloorY)*invRange;
			//beginVelX = 6.0*vMean*s*(1.0-s);
			beginVelX = Random(-vMean,vMean);
			beginVelY = Random(-vMean,vMean);
			ballA[nrBalls].initialize(nrBalls, CPair(beginPosX,beginPosY),
									  CPair(beginVelX,beginVelY),radius,mass);
			nrBalls++;
			if(nrBalls == maxNrBalls){stop = true;}
			beginPosX += maxRadius + 1.001*maxRadius;
			radius = Random(minRadius, maxRadius);
			//mass	  = (radius/maxRadius)*ballMass;
			mass = Random(minMass,maxMass);
		}
		beginPosY += maxRadius + 1.001*maxRadius;	
		beginPosX = inletX + 1.001*maxRadius;
	}

	beginPosX = tresholdX + 1.001*radius;
	beginPosY = outletFloorY + 1.001*radius;
	// the rest of the body
	/*while((beginPosY<ceilY-ballRadius)&&(!stop))
	{
		while((beginPosX<outletX-ballRadius)&&(!stop))
		{
			beginVelX = 0.1*vMean*Random();;
			beginVelY = 0.01*vMean*Random();
			ballA[nrBalls].initialize(nrBalls, CPair(beginPosX,beginPosY),
									  CPair(beginVelX,beginVelY));
			nrBalls++;
			if(nrBalls == maxNrBalls){stop = true;}
			beginPosX += ballRadius + 1.001*ballRadius;
		}
		beginPosY += ballRadius + 1.001*ballRadius;	
		beginPosX = tresholdX + 1.001*ballRadius;
	}*/
/*	radius = Random(minRadius, maxRadius);
	ballA[0].initialize(nrBalls,CPair(0.02,0.05),
						CPair(1,0), maxRadius);
	nrBalls++;
	radius = Random(minRadius, maxRadius);
	ballA[1].initialize(nrBalls,CPair(0.04,0.05),
						CPair(-1,0), minRadius);
	nrBalls++;
/*ballA[1].initialize(CPair(tresholdX+3*1.1*ballRadius,inletFloorY-1.1*ballRadius),
						CPair(0,0));
	nrBalls++;*/
	cout<<"Number of balls: "<<nrBalls<<endl;

	grid.update(ballA, nrBalls);
}

// draw the current configuration
// ------------------------------
void CContainer::draw()
{		
	// A // the grid in green
	/////// -----------------
	glColor3f(0,1,0);
	grid.draw();

	// B // the walls in blue
	/////// -----------------
	glColor3f(0,0,1);
	//DrawLine( inletX,	 ceilY,		   inletX,	  inletFloorY	);
	//DrawLine( inletX,	 inletFloorY,  tresholdX, inletFloorY	);
	//DrawLine( tresholdX, inletFloorY,  tresholdX, outletFloorY	);
	DrawLine( inletX,	 ceilY,		   inletX,	  outletFloorY	);
	DrawLine( inletX,	 outletFloorY, outletX,	  outletFloorY	);
	DrawLine( tresholdX, outletFloorY, outletX,	  outletFloorY	);
	DrawLine( outletX,	 outletFloorY, outletX,	  ceilY			);
	DrawLine( outletX,	 ceilY,		   inletX,	  ceilY			);

	// C // the balls in blue
	/////// -----------------
/*	glColor3f(0,0,1);
	for(int i=0; i<nrBalls; ++i)
	{
		ballA[i].draw();
	}
*/
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
//	checkBoundaries();
//	grid.imposeVelX(tresholdX-4.1*maxRadius,tresholdX);
//	grid.imposeVelX(inletX+30*maxRadius,inletX+40*maxRadius);
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
		CPair swarmVel(0,0);
		double IIswarmVelII;
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
						swarmVel = swarmVel + partnerVel;
						N = pPartner->getPos() - ballA[n].getPos(); // points from n to j
						distance = VectorLength(N);
						overlap = -(distance - ballA[n].getRadius() - pPartner->getRadius());
						if(overlap>=0) { 
			//				glColor3f(0,1,0); DrawPoint(ballA[n].getPos()+0.5*N);
							N = N/distance; // unit length now
							relVel = partnerVel - ballA[n].getVel();
							F = ballA[n].calculateF(k, c, overlap, N, relVel);
							ballA[n].addF(F); pPartner->addF(-F);
						//	ballA[n].bounce(pPartner);
						}
					}
				}
			}
		}
		IIswarmVelII = VectorLength(swarmVel);
		if(IIswarmVelII!=0){ swarmVel = swarmVel/IIswarmVelII;}
		else{swarmVel = CPair(0.70710678,0.70710678);}
		ballA[n].redirect(swarmVel);
	}
	
	// B // contact with a wall
	/////// -------------------
	for(i=0; i<nrBalls; ++i)
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
				//if(x>inletX+3*maxRadius)
			//	{ballA[i].setVelX(friction*ballA[i].getVel().x);}
			//	ballA[i].bounce();
			}
		}

/*		// contact with lower left wall
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
*/
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
		
				//if(x>inletX+3*maxRadius)
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
//	grid.checkBoundaries(ballA, nrBalls);
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