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
extern int activeTriangle;

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
	void initialize(int newNr, int newNrCorners, Corner* cornerA, int nrSides, int nrSlices);
	void	draw();
	void	drawFull();	
	void	calculateForce();
	double	calculateSurface();
	double	calculateVolume();
	double	getVolume(){ return volumeT;}
	double	getSurface(){ return surfaceT; }
	double	getPressure(){ return pressure; }
	void	save( ofstream outFile );
	int		getNr(){ return nr; }
	void	updateCenters();
	void	updateNormals();
private:
	int		nr;
	double	k;				// the proportional constant to regulate a constant surface
	double	c;				// the damping constant to regulate a constant surface
	Corner** pCornerA;
	CTriple	cellCenter;
	Spring** pSpringA;
	CTripleInt* triangleA;
	CTriple* normalA;
	int		nrCorners;
	int		nrSprings;
	int		nrSpring2s;
	int		nrTriangles;
	int		nrNormals;
	double	restVolume;
	double	volumeTminusOne;
	double	volumeT;
	double	restSurface;		// the surface at the beginning wich is the equilibrium (target)
	double	surfaceTminusOne;	// the surface at the previous timestep
	double	surfaceT;			// the surface at the current timestep
	double	thickness;			// thickness of the tissue / cell (120 µm for onion)
	double	pressure;			// the current pressure in the cell (N/m²)
	CTriple	calculateNormal(int nrVertices, CTriple* vertexA);
	CTriple calculateNormal(CTriple A, CTriple B, CTriple C);
	double	calculateSurfaceTriangle(CTriple A,CTriple B,CTriple C);
	void	drawPressureForce(CTriple A, CTriple B, CTriple C, CTriple Normal, double pressure);
	double	calculateVolumeTetrahedron(int triangleNr, CTriple A,CTriple B,CTriple C);
};

CTriple Cell::calculateNormal(CTriple A, CTriple B, CTriple C)
{ // A,B and C are the cornerpoints of a triangle, walking counterclockwise
	CTriple first, second;
	first = B - A;
	second = C - B;
	return NormVector(first,second);
}

CTriple	Cell::calculateNormal(int nrVertices, CTriple* vertexA)
{ // calculate the normal of a n-gon by the Newell method (from Hill, 2001)
	CTriple normal = ZEROTRIPLE;
	for(int i=0; i<nrVertices-1; ++i)
	{
		normal.x += (vertexA[i].y-vertexA[i+1].y)*(vertexA[i].z+vertexA[i+1].z);
		normal.y += (vertexA[i].z-vertexA[i+1].z)*(vertexA[i].x+vertexA[i+1].x);
		normal.z += (vertexA[i].x-vertexA[i+1].x)*(vertexA[i].y+vertexA[i+1].y);
	}
	normal.x += (vertexA[nrVertices-1].y-vertexA[0].y)*(vertexA[nrVertices-1].z+vertexA[0].z);
	normal.y += (vertexA[nrVertices-1].z-vertexA[0].z)*(vertexA[nrVertices-1].x+vertexA[0].x);
	normal.z += (vertexA[nrVertices-1].x-vertexA[0].x)*(vertexA[nrVertices-1].y+vertexA[0].y);
	normal = normal/VectorLength(normal);
	return normal;
}

void Cell::updateNormals()
{
	CTriple A,B,C;
	int i;
	for( i=0; i<nrNormals; ++i )
	{
		A = pCornerA[triangleA[i].x]->getPos();
		B = pCornerA[triangleA[i].y]->getPos();
		C = pCornerA[triangleA[i].z]->getPos();
		normalA[i] = CTriple(this->calculateNormal(A,B,C));
	}
}

// set the initial values
// ----------------------
void Cell::initialize(int newNr, int newNrCorners, Corner* cornerA, int nrSides, int nrSlices)
{
	nr			= newNr;
	// put pointers to the corners in place
	nrCorners	= newNrCorners;
	pCornerA	= new Corner*[nrCorners];
	nrTriangles	= 2*nrSides+2*(nrSlices-1)*nrSides;
	triangleA	= new CTripleInt[nrTriangles];
	nrNormals	= nrTriangles;
	normalA		= new CTriple[nrNormals];
	int i,j;
	// this is a simple copy since there is only one cell
	for(i=0; i<nrCorners; ++i)
	{
		pCornerA[i] = &cornerA[i];
	}
	this->updateCenters();
	// set the triangles, store the local numbers of the corners
	// get the faces the centerBottom with the cornerpoints of the bottom slice
	// faces are always defined in counterclockwise direction looking from outside
	int faceNr = 0;
	for(i=0; i<nrSides-1; ++i)
	{
		triangleA[faceNr] = CTripleInt(0,i+2,i+1);
		faceNr++;
	}
	triangleA[faceNr] = CTripleInt(0,1,nrSides);
	faceNr++;

	// calculate the numbers of the cornerpoints of the faces on the sides
	for(j=0; j<nrSlices-1; ++j)
	{
		for(i=0; i<nrSides-1; ++i)
		{
			triangleA[faceNr] = CTripleInt(1+i+j*nrSides,2+i+j*nrSides,2+i+(j+1)*nrSides);
			faceNr++;
			triangleA[faceNr] = CTripleInt(1+i+(j+1)*nrSides,1+i+j*nrSides,2+i+(j+1)*nrSides);
			faceNr++;
		}
		triangleA[faceNr] = CTripleInt(1+(nrSides-1)+j*nrSides,1+j*nrSides,1+(j+1)*nrSides);
		faceNr++;
		triangleA[faceNr] = CTripleInt(1+(nrSides-1)+(j+1)*nrSides,1+(nrSides-1)+j*nrSides,1+(j+1)*nrSides);
		faceNr++;
	}

	// add the top
	for(i=0; i<nrSides-1; ++i)
	{
		triangleA[faceNr] = CTripleInt(nrCorners-1,(nrCorners-1)-(i+2),(nrCorners-1)-(i+1));
		faceNr++;
	}
	triangleA[faceNr] = CTripleInt(nrCorners-1,nrCorners-1-1,nrCorners-1-nrSides);
	faceNr++;

	this->updateNormals();

	surfaceT = surfaceTminusOne = restSurface = this->calculateSurface();
	volumeT = volumeTminusOne = restVolume = this->calculateVolume();

//	thickness = 120e-6; // for unicellular onion tissue
	pressure = 0;

	{ // read from parameterFile.txt
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);

		iF.setScope("-parameters-");
		k	= iF.getDouble("kpressure");
		c	= 2*sqrt(k*cornerA[0].getMass()); 
		cout<<"Dpressure = "<<c<<endl;
	}

}

double Cell::calculateVolumeTetrahedron(int triangleNr, CTriple A,CTriple B,CTriple C)
{
	CTriple Normal = normalA[triangleNr];

	// we need the distance from the cellCenter to the plane defined by the triangle(A,B,C)
	CTriple D = A-cellCenter;
	double distance = D*Normal; // since N has unit length, the dot product is equal to the required distance now

	// calculate the area of the triangle
	double area = 0.5*VectorLength(VectorProd(B-A,C-B));

	// now we can calculate the volume of the tetrahedron(A,B,C,cellCenter)
	double volume = area*distance/3.0;
	return volume;
}

double Cell::calculateSurfaceTriangle(CTriple A,CTriple B,CTriple C)
{
	double area = 0.5*VectorLength(VectorProd(B-A,C-B));
	return area;
}

void Cell::drawPressureForce(CTriple A, CTriple B, CTriple C, CTriple Normal, double pressure)
{
	CTriple triangleCenter;
	triangleCenter = (A + B + C)/3.0;
	if(pressure>0)
	{	glColor3f(0,1,0);
	}else{
		glColor3f(1,0,0);
	}
		drawLine(triangleCenter ,triangleCenter+scaleFactorForForces*pressure*Normal);

}

double Cell::calculateVolume()
{
	double volume = 0;
	for(int i=0; i<nrTriangles; ++i)
	{
		CTriple A,B,C;
		A = pCornerA[triangleA[i].x]->getPos();
		B = pCornerA[triangleA[i].y]->getPos();
		C = pCornerA[triangleA[i].z]->getPos();

		volume += this->calculateVolumeTetrahedron(i,A,B,C);
	}
	return volume;
}

// calculate the force due to the internal pressure, the spring forces are calculated via container
// -------------------------------------------
void Cell::calculateForce()
{
	int i;
	double surface;
	CTriple Normal, F;

	// forces exerted by the internal pressure
	// ---------------------------------------
	surfaceTminusOne = surfaceT;
	surfaceT = this->calculateSurface();
	volumeTminusOne = volumeT;
	volumeT = this->calculateVolume();

	// 'surface velocity'/dt to make it [m²/s]
	pressure = -(k*(volumeT-restVolume)+c*(volumeT-volumeTminusOne)/dt);

//	glColor4f(1,0,0,1);
//	drawDisc(cellCenter,0.00008);

	for(i=0; i<nrTriangles; ++i)
	{
		CTriple A,B,C;
		A = pCornerA[triangleA[i].x]->getPos();
		B = pCornerA[triangleA[i].y]->getPos();
		C = pCornerA[triangleA[i].z]->getPos();
		
		Normal = normalA[i];
		surface = this->calculateSurfaceTriangle(A,B,C);

		F = pressure*surface*Normal;

		if(showForces){
			this->drawPressureForce(A,B,C,Normal,pressure);
		}
		// divide the pressure force over the three corners
		pCornerA[triangleA[i].x]->addFtoAccP(F/3.0);
		pCornerA[triangleA[i].y]->addFtoAccP(F/3.0);
		pCornerA[triangleA[i].z]->addFtoAccP(F/3.0);
	}
}



void Cell::draw() 
{ 
	int i;
	CTriple pos;
	glBegin(GL_TRIANGLES); 
	for( i=0; i<nrTriangles; ++i )
	{
//		glColor4f(i/(double)(nrTriangles),1-i/(double)(nrTriangles),0.5,0.3);
		glColor4f(0,1,0,0.3);
//		if(i==activeTriangle){glColor4f(0,0,1,1);}
		pos = pCornerA[triangleA[i].x]->getPos(); 
		glVertex3f(pos.x,pos.y,pos.z);
		pos = pCornerA[triangleA[i].y]->getPos(); 
		glVertex3f(pos.x,pos.y,pos.z);
		pos = pCornerA[triangleA[i].z]->getPos(); 
		glVertex3f(pos.x,pos.y,pos.z);
	}
	glEnd();

/*	CTriple begin, end;
	for( i=0; i<nrNormals; ++i )
	{
		if(i==activeTriangle){
			glColor4f(1,0,0,1);
			begin = pCornerA[triangleA[i].x]->getPos();
			end   = begin + scaleFactorForForces*normalA[i];
			drawLine(begin, end);
			begin = pCornerA[triangleA[i].y]->getPos();
			end   = begin + scaleFactorForForces*normalA[i];
			drawLine(begin, end);
			begin = pCornerA[triangleA[i].z]->getPos();
			end   = begin + scaleFactorForForces*normalA[i];
			drawLine(begin, end);
		}
	}
*/
}

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

void Cell::updateCenters()
{
	cellCenter = ZEROTRIPLE;
	int i;
	for(i=0; i<nrCorners; ++i){
		cellCenter = cellCenter + pCornerA[i]->getPos();
	}
	cellCenter = cellCenter/(double)(nrCorners);
}

// calculate the surface, supposing each next triangle contributes a positive part
// -------------------------------------------------------------------------------
double Cell::calculateSurface()
{
	double	surface	= 0;	// the total surface
	int i;
	CTriple	A, B, C;
	CTriple IvectorProdI;
	A = pCornerA[0]->getPos();

	for(i=0; i<nrCorners-1; ++i)
	{
		B = pCornerA[i]->getPos()-A;
		C = pCornerA[i+1]->getPos()-A;
		IvectorProdI = VectorProd(B,C);
		
		surface += 0.5*Dir(IvectorProdI)*VectorLength(IvectorProdI);
	}

	return surface; 
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
}




//////////////////////////////////////////////////////////////////
#endif															//
// end of Cell.h												//
//////////////////////////////////////////////////////////////////