//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include <iostream> 

#include "General.h"
#include "Corner.h"
#include "Cell.h"
#include "InputFile.h"
#include "BoundaryConditions.h"

extern bool	  showForces;	// from main.cpp
extern int	  forceMode;	// from main.cpp

/**
 * The Container class. It contains an array of:
 * - all the corners involved in the simulation
 * - the cells
 * The Container should be declared as a global object in main() 
 * to control the simulation
 */
class Container
{	
public:	
	Container();
	~Container(){};
	void	draw();
	void	move( double dt );
	void	calculateForces();
	void	applyBC(int mode){ BC.apply(mode); }
	void	save();
	void	saveData( double time );
private:
	double	calculateSurface();
	Cell   *cellA;
	Corner *cornerA;
	Spring *springA;
	Spring *secSpringA;
	BoundaryConditions BC;
	int		nrCells;
	int		nrCorners;
	int		nrSprings;
	int		nrSecSprings;
	double	initialLength;	// initial length of the whole configuration in m
	double	length;			// the current length in m
	double	strain;			// initially there is no strain of the tissue
	double	stress;			// the stress measured at the north face
	double	initialSurface;	// in m²
	double	initialCellSurface; // in m²
	double	surface;			// the current surface in m²
	double	pressure;		// the average pressure in a cell (Pa)
	double	pressureMin;
	double	pressureMax;
};



// set the beginconditions
// -----------------------
Container::Container( )
{
	const char* nameOut = "data/forceDeformation.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	outFile.close();

	double newCornerMass;
	double newK, newC, newSecK, newSecC;
	{
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);
		newCornerMass	= iF.getDouble("mass"); 
		newK			= iF.getDouble("k");
		newC			= 2*sqrt(newK*newCornerMass);
		newSecK			= iF.getDouble("ksecundary");
		newSecC			= 2*sqrt(newSecK*newCornerMass);
		cout<<"DspringCell = "<<newC<<endl;
		cout<<"DspringPlane = "<<newSecC<<endl;
	}
	

	int i;
	const char* nameIn = "data/inputFile.txt";
	InputFile iF(nameIn);

	// read the coordinates of the cornerpoints
	iF.setScope("-cornerpoints-");
	nrCorners = iF.getInt("nrcorners");
	cornerA = new Corner[nrCorners];
	CTriple	P;
	for(i=0; i<nrCorners; ++i)
	{
		int j = iF.getInt();
		if(j!=i+1){cout<<"error in inputFile"<<endl;}
		P.x		= iF.getDouble();
		P.y		= iF.getDouble();
		P.z		= iF.getDouble();

		cornerA[i].initialize(j,P,newCornerMass);
	}

	// read the springs
	iF.setScope("-cellboundarysprings-");
	nrSprings = iF.getInt("nrsprings");
	springA = new Spring[nrSprings];
	Corner *begin, *end;
	for(i=0; i<nrSprings; ++i)
	{
		int j = iF.getInt();
		if(j!=i+1){cout<<"error in inputFile"<<endl;}
		begin	= &cornerA[iF.getInt()-1];
		end		= &cornerA[iF.getInt()-1];
		springA[i].initialize(j,newK,newC,begin,end);
	}

	// read the secundary springs
	iF.setScope("-planeboundarysprings-");
	nrSecSprings = iF.getInt("nrsprings");
	secSpringA = new Spring[nrSecSprings];
	for(i=0; i<nrSecSprings; ++i)
	{
		int j = iF.getInt();
		if(j!=i+1){cout<<"error in inputFile"<<endl;}
		begin	= &cornerA[iF.getInt()-1];
		end		= &cornerA[iF.getInt()-1];
		secSpringA[i].initialize(j,newSecK,newSecC,begin,end);
	}


	// read the configuration
	iF.setScope("-configuration-");
	nrCells = iF.getInt("nrcells");
	int  *nrCornersA	= new int [nrCells];
	int **cornerNrsA	= new int*[nrCells];
	int  *nrSpringsA	= new int [nrCells];
	int **springNrsA	= new int*[nrCells];
	for(i=0; i<nrCells; ++i)
	{
		int nrCorners = iF.getInt("nrcorners");
		nrCornersA[i] = nrCorners;
		cornerNrsA[i] = new int[nrCorners];
		cornerNrsA[i][0] = iF.getInt("numbers");
		int j;
		for(j=1; j<nrCorners; ++j)
		{
			cornerNrsA[i][j] = iF.getInt();
		}
		iF.setScopeBegin(iF.location);
		int nrSprings = iF.getInt("nrsprings");
		nrSpringsA[i] = nrSprings;
		springNrsA[i] = new int[nrSprings];
		springNrsA[i][0] = iF.getInt("numbers");
		for(j=1; j<nrSprings; ++j)
		{
			springNrsA[i][j] = iF.getInt();
		}
		iF.setScopeBegin(iF.location);
	}

	// initialize all the cells with their respective corners
	cellA = new Cell[nrCells];
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].initialize(i+1, nrCornersA[i], cornerNrsA[i], cornerA,
							     nrSpringsA[i], springNrsA[i], springA);
	}
	delete[] nrCornersA;
	delete[] cornerNrsA;
	delete[] nrSpringsA;
	delete[] springNrsA;

	//set the window such that the last corner fits in
	// further adjusted in BoundaryConditions::initialize()
	extern double xmin,xmax,ymin,ymax;
	CTriple max = cornerA[nrCorners-1].getPos();
	CTriple min = cornerA[0].getPos();
	CTriple pos = max - min;
	if(pos.x<pos.y){
		xmin = min.x-0.5*(pos.y-pos.x);
		xmax = min.x+pos.x+0.5*(pos.y-pos.x);
		ymin = min.y;
		ymax = min.y+pos.y;	
	} else {
		xmin = min.x;
		xmax = min.x+pos.x;
		ymin = min.y-0.5*(pos.x-pos.y);
		ymax = min.y+pos.y+0.5*(pos.x-pos.y);
	}

	// initialize the BC
	BC.initialize(cornerA, springA, secSpringA);
	initialLength = BC.getLength();
	// initialize the measured quantities
	strain = 0;
	stress = 0;
	initialSurface = surface = this->calculateSurface();
	initialCellSurface = cellA[0].calculateSurface();

}

// draw the current configuration
// ------------------------------
void Container::draw()
{
	int i=0;
	glColor3f(0,0,1);
	for( i=0; i<nrCorners; ++i ){
		drawDisc(cornerA[i].getPos(), 0.00002);
	}
	glColor3f(0,0.8,1);
	for( i=0; i<nrSprings; ++i ){ springA[i].draw(); }
	glColor3f(0,0,1);
	for( i=0; i<nrSecSprings; ++i ){ secSpringA[i].draw(); }
/**/
}

// move the corners by integrating the equation of motion
// ------------------------------------------------------
void Container::move(double dt)
/*
 * The corners updated their accelerations in the calculateForces step
 * now their positions will be updated
 */
{
	// make sure the north and south face do not move
	int i=0;
	// integrate the accelerations of the corners
	for(i=0; i<nrCorners; ++i)
	{
		cornerA[i].move(dt);
	}
}


// calculate all the forces acting on the corners
// ----------------------------------------------
void Container::calculateForces()
{
	int i;
	// the spring forces
	for(i=0; i<nrSprings; ++i)
	{
		springA[i].calculateForce();
	}

	// the secundary spring forces
	for(i=0; i<nrSecSprings; ++i)
	{
		secSpringA[i].calculateForce();
	}
	
	// the pressure forces
	if(forceMode == 1) {
		for(i=0; i<nrCells; ++i)
		{
			cellA[i].calculateForce();
		}
	}
}

// calculate the total surface of the tissue
// -----------------------------------------
double Container::calculateSurface()
{
	double surface = 0;
	int i;
	for(i=0; i<nrCells; ++i)
	{
		surface += cellA[i].calculateSurface(); // use getSurface when surf is kept constant
	}
	return surface;
}

void Container::saveData(double time)
{
	// update the quantities
	BC.setStrain(100.0*((cornerA[15].getPos().x-cornerA[13].getPos().x)-0.00012)/0.00012);
	strain	= BC.getStrain();
	stress	= BC.getStress();
	surface = this->calculateSurface();
	double pressureSum = 0;
	double pressureMin = 0;
	double pressureMax = 0;

	for(int i=0; i<nrCells; ++i)
	{
		pressure = cellA[i].getPressure();
		if(pressure<pressureMin){pressureMin = pressure;}
		if(pressure>pressureMax){pressureMax = pressure;}
		pressureSum += pressure;
	}
	pressure = pressureSum/(double)(nrCells);

	const char* nameOut = "data/forceDeformation.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::app);
	if (! outFile){}
	outFile	<<time<<" "
			<<strain<<" "
			<<stress<<" "
			<<surface<<" "
			//<<100.0*(surface-initialSurface)/initialSurface<<" "
			<<BC.getWidth()<<" "
			<<pressure<<" "
			<<pressureMin<<" "
			<<pressureMax<<" ";

			for(i=0; i<nrCells; ++i)
			{		
				outFile<<100.0*(cellA[i].calculateSurface()-initialCellSurface)/initialCellSurface<<" ";
			}
	outFile<<endl;
	outFile.close();
//	cout<<"strain = "<<strain<<endl;
}

// save the container and everything else to be able to resume calculation from these data
// ---------------------------------------------------------------------------------------
void Container::save()
{
	const char* nameOut = "data/savedInputFile.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	if (! outFile){}

	// save the cornerpoints
	outFile<<"-cornerpoints-"<<endl;
	outFile<<"nrcorners "<<nrCorners<<endl;
	int i;
	for(i=0; i<nrCorners; ++i)
	{
		cornerA[i].save(outFile);
	}
	outFile<<"-end-cornerpoints-"<<endl;
	outFile<<endl;
	
	// save the springs
	outFile<<"-springs-"<<endl;
	outFile<<"nrsprings "<<nrSprings<<endl;
	for(i=0; i<nrSprings; ++i)
	{
		springA[i].save(outFile);
	}
	outFile<<"-end-springs-"<<endl;
	outFile<<endl;

	// save the configuration
	outFile<<"-configuration-"<<endl;
	outFile<<"nrcells "<<nrCells<<endl;
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].save(outFile);	
	}
	outFile<<"-end-configuration-"<<endl;
	outFile<<endl;

	// save the boundaryconditions
	BC.save(outFile);
	outFile<<endl;
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////