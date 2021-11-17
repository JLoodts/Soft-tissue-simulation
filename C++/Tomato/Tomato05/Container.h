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
	void	applyBC(int mode){ BC.apply(mode,nrCorners,cornerA); }
	void	save();
	void	saveData( double time );
private:
	double	calculateSurface();
	double	calculateVolume();
	Cell   *cellA;
	Corner *cornerA;
	Spring *springA;
	BoundaryConditions BC;
	int		nrCells;
	int		nrCorners;
	int		nrSprings;
	double	initialHeight;	// initial height of the whole configuration (it est cell) in m
	double	length;			// the current length in m
	double	strain;			// initially there is no strain of the tissue
	double	forceY;			// the force.y measured from the southern hemisphere
	double	initialSurface;	// in m²
	double	initialCellSurface; // in m²
	double	volume;
	double	initialVolume;
	double	initialCellVolume;
	double	surface;			// the current surface in m²
	double	pressure;		// the average pressure in a cell (Pa)
	double	pressureMin;
	double	pressureMax;
	void	initialize(double radius, CTriple center, int nrSlices, int nrSides);
};

// draw the current configuration
// ------------------------------
void Container::draw()
{
	int i=0;
/*	glColor3f(0,0,1);
	for( i=0; i<nrCorners; ++i ){	
		glColor3f(i/(double)(nrCorners),0,1-i/(double)(nrCorners));
//		drawDisc(cornerA[i].getPos(), 0.00002);
		cornerA[i].draw();
	}
*/	BC.draw();
/*	glColor3f(0,0.8,1);
	for( i=0; i<nrSprings; ++i ){ springA[i].draw(); }
*/	glColor3f(0,0,1);
	for( i=0; i<nrCells; ++i ){ cellA[i].draw(); }

}

// set the beginconditions
// -----------------------
Container::Container( )
{
	this->initialize(62e-6,CTriple(0,0,0),40,45);
}

void Container::initialize(double radius, CTriple center, int nrSlices, int nrSides)
{
	// nrSlices: the number of circular discs of cornerpoints created to form the sphere
	//			 the north and south point are not taken into account for this

	const char* nameOut = "data/forceDeformation.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	outFile.close();

	double newCornerMass;
	double newK, newC;
	{
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);
		newCornerMass	= iF.getDouble("mass"); 
		newK			= iF.getDouble("k");
		newC			= 2*sqrt(newK*newCornerMass);
		cout<<"DspringCell = "<<newC<<endl;
	}
	

	int i,j;

	// calculate the coordinates of the cornerpoints
	nrCorners = nrSides*(nrSlices)+2;
	cornerA = new Corner[nrCorners];
	CTriple	P, centerBottom, centerJ;
	centerBottom = center - CTriple(0,radius,0);
	double ySlice;
	double alfa = 2.0*PI/(double)(nrSides);
	double beta = PI/(double)(nrSlices+1);
	double radiusJ;
	int cornerNr = 0;
	cornerA[cornerNr].initialize(cornerNr+1,centerBottom,newCornerMass);
	cornerNr++;
	for(j=1; j<=nrSlices; ++j)
	{
		double corner = (double)(j)*beta;
		radiusJ = radius*sin(corner);
		ySlice	= radius*(1-cos(corner));

		centerJ = centerBottom + CTriple(0,ySlice,0);
		for(i=0; i<nrSides; ++i)
		{
			P.x		= centerJ.x + radiusJ*cos((double)(i)*alfa);
			P.y		= centerJ.y;
			P.z		= centerJ.z - radiusJ*sin((double)(i)*alfa);

			cornerA[cornerNr].initialize(cornerNr+1,P,newCornerMass);
			cornerNr++;
		}
	}
	cornerA[cornerNr].initialize(cornerNr+1,centerBottom+CTriple(0,2.0*radius,0),newCornerMass);
	cout<<cornerNr+1<<" corners were created"<<endl;
	// calculate the springs
	nrSprings = 2*nrSides + nrSlices*nrSides + (nrSlices-1)*3*nrSides; 
	springA = new Spring[nrSprings];
	Corner *begin, *end;
	int springNr = 0;
	// connect the centerBottom with the cornerpoints of the bottom slice
	for(i=0; i<nrSides; ++i)
	{
		begin	= &cornerA[0];
		end		= &cornerA[1+i];
		springA[springNr].initialize(springNr+1,newK,newC,begin,end);
		springNr++;
	}

	// connect all the corners contained in the same slice and this for each slice
	for(j=0; j<nrSlices; ++j)
	{
		for(i=0; i<nrSides-1; ++i)
		{
			begin	= &cornerA[1+i+j*nrSides];
			end		= &cornerA[1+i+j*nrSides+1];
			springA[springNr].initialize(springNr+1,newK,newC,begin,end);
			springNr++;
		}
		begin	= &cornerA[1+(nrSides-1)+j*nrSides];
		end		= &cornerA[1+(0)+j*nrSides];
		springA[springNr].initialize(springNr+1,newK,newC,begin,end);
		springNr++;
	}

	// connect all the corners of neighbouring slices (springs parallel to y-axis)
	for(j=0; j<nrSlices-1; ++j)
	{
		for(i=0; i<nrSides; ++i)
		{
			begin	= &cornerA[1+i+j*nrSides];
			end		= &cornerA[1+i+(j+1)*nrSides];
			springA[springNr].initialize(springNr+1,newK,newC,begin,end);
			springNr++;
		}
	}

	// connect all the corners of neighbouring slices (springs are diagonal of faces)
	for(j=0; j<nrSlices-1; ++j)
	{
		for(i=0; i<nrSides-1; ++i)
		{
			begin	= &cornerA[1+i+j*nrSides];
			end		= &cornerA[1+i+(j+1)*nrSides+1];
			springA[springNr].initialize(springNr+1,newK,newC,begin,end);
			springNr++;
			begin	= &cornerA[1+i+(j+1)*nrSides];
			end		= &cornerA[1+i+j*nrSides+1];
			springA[springNr].initialize(springNr+1,newK,newC,begin,end);
			springNr++;
		}
		begin	= &cornerA[1+(nrSides-1)+j*nrSides];
		end		= &cornerA[1+(0)+(j+1)*nrSides];
		springA[springNr].initialize(springNr+1,newK,newC,begin,end);
		springNr++;
		begin	= &cornerA[1+(nrSides-1)+(j+1)*nrSides];
		end		= &cornerA[1+(0)+j*nrSides];
		springA[springNr].initialize(springNr+1,newK,newC,begin,end);
		springNr++;
	}

	// connect the centerTop with the cornerpoints of the top slice
	for(i=0; i<nrSides; ++i)
	{
		begin	= &cornerA[nrCorners-1];
		end		= &cornerA[nrCorners-2-i];
		springA[springNr].initialize(springNr+1,newK,newC,begin,end);
		springNr++;
	}

	// initialize all the cells with their respective corners
	nrCells = 1;
	cellA = new Cell[nrCells];
	cellA[0].initialize(1, nrCorners, cornerA, nrSides, nrSlices);

	//set the window such that the last corner fits in
	// further adjusted in BoundaryConditions::initialize()
	extern double xmin,xmax,ymin,ymax;
	CPair max = CPair(center.x,center.y) + radius*CPair(1,1);
	CPair min = CPair(center.x,center.y) - radius*CPair(1,1);
	CPair pos = max - min;
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
	BC.initialize(nrCorners, cornerA, springA);
	initialHeight = BC.getHeight();

	// initialize the measured quantities
	strain = 0;
	forceY = 0;
	initialSurface = surface = this->calculateSurface();
	initialCellSurface = cellA[0].calculateSurface();
	initialVolume = volume = this->calculateVolume();
	initialCellVolume = cellA[0].calculateVolume();
	initialCellVolume = volume;
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
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].updateCenters();
		cellA[i].updateNormals();
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

double Container::calculateVolume()
{
	double volume = 0;
	int i;
	for(i=0; i<nrCells; ++i)
	{
		volume += cellA[i].calculateVolume();
	}
	return volume;
}

void Container::saveData(double time)
{
	// update the quantities
	BC.setStrain(100.0*((cornerA[15].getPos().x-cornerA[13].getPos().x)-0.00012)/0.00012);
	strain	= BC.getStrain();
	forceY	= BC.getForceY(nrCorners, cornerA);
	surface = this->calculateSurface();
	volume	= this->calculateVolume();
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
			<<forceY<<" "
			<<volume<<" "
			<<100.0*(volume-initialVolume)/initialVolume<<" "
			//<<100.0*(surface-initialSurface)/initialSurface<<" "
			//<<BC.getWidth()<<" "
			<<pressure<<" "
			<<pressureMin<<" "
			<<pressureMax<<" ";

			outFile<<surface<<" ";
			/*for(i=0; i<nrCells; ++i)
			{		
				outFile<<100.0*(cellA[i].calculateVolume()-initialCellVolume)/initialCellVolume<<" ";
			}*/
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
}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////