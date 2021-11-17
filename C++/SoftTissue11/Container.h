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
extern double loadMax;		// from main.cpp
extern int	  activeCell;	// from main.cpp
extern double pullForceX;	// from main.cpp
extern double pullForceY;	// from main.cpp
extern double displacementSpeed; // from main.cpp


class CContainer
/*
 * Contains an array of the corners which are linked to one another
 * each corner has its own dynamics which is determined by spring forces
 * connecting it to its two neighbours and a pressure force
 */
{	
public:	
	CContainer					();
	~CContainer					(){};
	void	continuePreviousCalculation();
	void	draw				();
	void	move				(double dt);
	void	calculateForces		(double time, double k, double c, double displacement, 
								 double& elongation, double& force);
	void	pullUp();
	double	getVolume(){ return cellA[activeCell].calculateSurface();}
	double	getNorthPos();
	double	getEastPos();
	double	getSouthPos();
	double	getWestPos();
	void	setNorthPos(double northPosY);
	void	setEastPos(double eastPosX);
	void	setSouthPos(double southPosY);
	void	setWestPos(double westPosX);
	void	setPullForceX();
	double	getDeformation(){return this->getNorthPos() - this->getSouthPos();}
	void	save		();
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
	int		nrMeasureStress;
	double	initialLength;
	double	thickness;
	CCorner *measureBeginCorner;
	CCorner *measureEndCorner;
	CCorner **northFace;
	CCorner **southFace;
	CCorner **eastFace;
	CCorner **westFace;
	CSpring **measureStressA;
};


void CContainer::continuePreviousCalculation()
{
	const char* nameIn = "data/inputFile.txt";
	InputFile iF(nameIn);
	// initialize all the cells with the information strored in the inputfile
	iF.setScope("-springsandpressure-");
	for(int i=0; i<nrCells; ++i)
	{
		int cellNr = iF.getInt("cellnr");
		cellA[cellNr-1].setInitialConditions(iF.getDouble("ratiopxv"),iF.getDouble("restpressure"));
		cellA[cellNr-1].setRestLength(0, iF.getDouble("restlength"));
		for(int i=1; i<cellA[cellNr-1].getNrCorners(); ++i)
		{
			cellA[cellNr-1].setRestLength(i, iF.getDouble());
		}
		iF.setScopeBegin();
	}
}

// set the beginconditions
// -----------------------
CContainer::CContainer( )
{
	const char* nameOut = "data/forceDeformation.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	outFile.close();

	double newCornerMass;
	double newPressure;
	{
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);
		newCornerMass = iF.getDouble("mass"); 
		newPressure = iF.getDouble("pressure");
	}
	

	int i;
	const char* nameIn = "data/inputFile.txt";
	InputFile iF(nameIn);
	// read the coordinates of the cornerpoints
	iF.setScope("-cornerpoints-");
	nrCorners = iF.getInt("nrcorners");
	cornerA = new CCorner[nrCorners];
	CPair	P;
	double scalingFactor = 1.0;//0.00428;
	for(i=0; i<nrCorners; ++i)
	{
		int j = iF.getInt();
		if(j!=i+1){cout<<"error in inputFile"<<endl;}
		P.x		= scalingFactor*iF.getDouble();
		P.y		= scalingFactor*iF.getDouble();
		cornerA[i].initialize(i+1,P,newCornerMass);
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
		cellA[i].initialize(i+1, nrCornersA[i], cornerNrsA[i], cornerA, newPressure);
	}
	delete[] nrCornersA;
	delete[] cornerNrsA;
	

	// set the north and south positions correctly
//	initialLength = this->getNorthPos()-this->getSouthPos();
	initialLength = this->getEastPos()-this->getWestPos();

	// read the walls in which the stress has to be measured
	iF.setScope("-measurestress-");
	int cellNr, wallNr, cornerNr;
	cornerNr = iF.getInt("beginmeasurecorner");
	measureBeginCorner = &cornerA[cornerNr-1];
	cornerNr = iF.getInt("endmeasurecorner");
	measureEndCorner = &cornerA[cornerNr-1];
	thickness = iF.getDouble("thickness");
	
	nrMeasureStress = iF.getInt("nrwalls");
	measureStressA = new CSpring*[nrMeasureStress];
	for(i=0; i<nrMeasureStress; ++i)
	{
		cellNr = iF.getInt("cellnr");
		wallNr = iF.getInt("wallnr");
		measureStressA[i] = cellA[cellNr-1].getWallP(wallNr-1);
		iF.setScopeBegin(iF.location);
	}
}

void CContainer::save()
{
	const char* nameOut = "data/savedInputFile.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	if (! outFile){}
	outFile<<"-cornerpoints-"<<endl;
	outFile<<"nrcorners "<<nrCorners<<endl;
	for(int i=0; i<nrCorners; ++i)
	{

		outFile<<cornerA[i].getNr()<<" "<<cornerA[i].getPos().x<<" "<<cornerA[i].getPos().y<<endl;
	}
	outFile<<"-end-cornerpoints-"<<endl;
	outFile<<endl;
	outFile<<"-configuration-"<<endl;
	outFile<<"nrcells "<<nrCells<<endl;
	for(i=0; i<nrCells; ++i)
	{

		outFile<<cellA[i].getNr()<<" "<<" nrcorners "<<cellA[i].getNrCorners()<<" numbers ";
		for(int j=0; j<cellA[i].getNrCorners(); ++j)
		{
			outFile<<cellA[i].getCornerNr(j)<<" ";
		}
		outFile<<endl;			
	}
	outFile<<"-end-configuration-"<<endl;
	outFile<<endl;
	outFile<<"-northface-"<<endl;
	outFile<<"nrcorners "<<nrNorthFace<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrNorthFace; ++i)
	{
		outFile<<northFace[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-northface-"<<endl;
	outFile<<endl;

	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrEastFace<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrEastFace; ++i)
	{
		outFile<<eastFace[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-eastface-"<<endl;
	outFile<<endl;

	outFile<<"-southface-"<<endl;
	outFile<<"nrcorners "<<nrSouthFace<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrSouthFace; ++i)
	{
		outFile<<southFace[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-southface-"<<endl;
	outFile<<endl;

	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrWestFace<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrWestFace; ++i)
	{
		outFile<<westFace[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-westface-"<<endl;
	outFile<<endl;

	// save the restlength and pressure data
	outFile<<"-springsandpressure-"<<endl;
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].save(outFile);
		outFile<<endl;			
	}
	outFile<<"-end-springsandpressure-"<<endl;
	outFile.close();
}



// draw the current configuration
// ------------------------------
void CContainer::draw()
{	//draw the floor
	//glColor3f(0,0,1);
	//drawLine(CPair(0,floorHeight),CPair(1,floorHeight));
	int i=0;
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].draw();
	}
/*	for(i=0; i<nrCorners; ++i)
	{
		cornerA[i].draw();
	}
	glColor3f(1,0,0);
/*	extern int activeCell;
	cellA[activeCell].drawFull();
*/
}

// move the corners by integrating the equation of motion
// ------------------------------------------------------
void CContainer::move(double dt)
/*
 * The corners updated their accelerations in the calculateForces step
 * now their positions will be updated
 */
{
	// make sure the north and south face do not move
	int i=0;
	loadMax = 1e-6;
/*	for(i=0; i<nrNorthFace; ++i)
	{
		northFace[i]->setVel(ZEROPAIR);
		northFace[i]->setAcc(ZEROPAIR);
		northFace[i]->setVolAcc(ZEROPAIR);
	}
	for(i=0; i<nrSouthFace; ++i)
	{
		southFace[i]->setVel(ZEROPAIR);
		southFace[i]->setAcc(ZEROPAIR);
		southFace[i]->setVolAcc(ZEROPAIR);
	}

*/	

	for(i=0; i<nrEastFace; ++i)
	{
		eastFace[i]->setVel(ZEROPAIR);
		eastFace[i]->setAcc(ZEROPAIR);
		eastFace[i]->setVolAcc(ZEROPAIR);
	}
	for(i=0; i<nrWestFace; ++i)
	{
		westFace[i]->setVel(ZEROPAIR);
		westFace[i]->setAcc(ZEROPAIR);
		westFace[i]->setVolAcc(ZEROPAIR);
	}
	// integrate the accelerations of the corners
	for(i=0; i<nrCorners; ++i)
	{
		cornerA[i].move(dt);
	}
}

double	CContainer::getNorthPos()
{
	int i;
	double northPosY = 0;
	for(i=0; i<nrNorthFace; ++i)
	{
		northPosY += northFace[i]->getPos().y;
	}
	return northPosY = northPosY/(double)(nrNorthFace);
}

double	CContainer::getEastPos()
{
	int i;
	double eastPosX = 0;
	for(i=0; i<nrEastFace; ++i)
	{
		eastPosX += eastFace[i]->getPos().x;
	}
	return eastPosX = eastPosX/(double)(nrEastFace);
}

double	CContainer::getSouthPos()
{
	int i;
	double southPosY = 0;
	for(i=0; i<nrSouthFace; ++i)
	{
		southPosY += southFace[i]->getPos().y;
	}
	return southPosY = southPosY/(double)(nrSouthFace);
}

double	CContainer::getWestPos()
{
	int i;
	double westPosX = 0;
	for(i=0; i<nrWestFace; ++i)
	{
		westPosX += westFace[i]->getPos().x;
	}
	return westPosX = westPosX/(double)(nrWestFace);
}

void	CContainer::setNorthPos(double northPosY)
{
	int i;
	for(i=0; i<nrNorthFace; ++i)
	{
		northFace[i]->setPosY(northPosY);
	}
}

void	CContainer::setEastPos(double eastPosX)
{
	int i;
	for(i=0; i<nrEastFace; ++i)
	{
		eastFace[i]->setPosX(eastPosX);
	}
}

void	CContainer::setSouthPos(double southPosY)
{
	int i;
	for(i=0; i<nrSouthFace; ++i)
	{
		southFace[i]->setPosY(southPosY);
	}
}

void	CContainer::setWestPos(double westPosX)
{
	int i;
	for(i=0; i<nrWestFace; ++i)
	{
		westFace[i]->setPosX(westPosX);
	}

}

void CContainer::setPullForceX()
{
	int i;
	for(i=0; i<nrEastFace; ++i)
	{
		eastFace[i]->addF(CPair(pullForceX,0));
		//eastFace[i]->setVolAcc(ZEROPAIR);
		//eastFace[i]->setVel(ZEROPAIR);
	}	
	for(i=0; i<nrWestFace; ++i)
	{
		westFace[i]->addF(CPair(-pullForceX,0));
		//westFace[i]->setVolAcc(ZEROPAIR);
		//westFace[i]->setVel(ZEROPAIR);
	}
}

// calculate all the forces acting on the corners
// ----------------------------------------------
void CContainer::calculateForces(double time, double k, double c, double displacement,
								   double& elongation, double& force)
/**
 * Calculate the force and modify the accelerations of the corners
 * N points from the begincorner to the endcorner (is CCW) and the
 * resulting F acts on the respective corners
 */
{
	int i;

	// set the positions of the corners in north and south  correctly
//	this->setNorthPos(this->getNorthPos()+0.5*displacement);
//	this->setSouthPos(this->getSouthPos()-0.5*displacement);
	this->setEastPos(this->getEastPos()+0.5*displacement);
	for(i=0; i<nrEastFace; ++i)
	{
		eastFace[i]->setVel(CPair(displacementSpeed,0));
	}
	this->setWestPos(this->getWestPos()-0.5*displacement);
	for(i=0; i<nrWestFace; ++i)
	{
		westFace[i]->setVel(CPair(-displacementSpeed,0));
	}
//	this->setPullForceX();

	// A // the corners are connected by springs
	/////// ------------------------------------
	for(i=0; i<nrCells; ++i)
	{
//		if((i==3)||(i==4)) {mode=1;// also pressure and secondary springs
//		} else {mode=2;} // only the springs of the walls
		// forceMode = 1 pressure + spring
		//			!= 1 spring
		cellA[i].calculateF(k,c);
	}


	// calculate the force acting in the walls of a cut through
	double F = 0;
	for(i=0; i<nrMeasureStress; ++i)
	{
		F = F + measureStressA[i]->getLoad();
		glColor3f(1,0,0);
		measureStressA[i]->draw();
	}
	double width = Distance(measureBeginCorner->getPos(),measureEndCorner->getPos());
	double area = thickness*width;

	elongation = 100*(this->getEastPos()-this->getWestPos()-initialLength)/(double)initialLength;
	force = F/area; // divide by the area in m²

/*	CPair F;
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

	// get loadMax right for this timestep
	
	for(i=0; i<nrCells; ++i)
	{
		//if(cornerA[i].getLoad()>loadMax){ loadMax = cornerA[i].getLoad();}
		if(cellA[i].getLoadMax()>loadMax){ loadMax = cellA[i].getLoadMax();}
	}
*/
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////