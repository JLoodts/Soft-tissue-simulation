//////////////////////////////////////////////////////////////////
#ifndef BOUNDARYCONDITIONS_H									//
#define BOUNDARYCONDITIONS_H									//
// begin of BoundaryConditions.h								//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include "InputFile.h"

extern double dt; // from main.cpp
extern int	  timeSteps; // from main.cpp

/**
 *  The BoundaryConditions class. 
 *	
 *	Handles the BC that are applied
 *	Wether the tissue is stretched or relaxed and at which velocity this happens.
 *	Has pointers to the corners in the east and the west
 */
class BoundaryConditions
{	
public:
	BoundaryConditions(){}
	~BoundaryConditions(){}
	void		initialize(int nrCorners, Corner* cornerA, Spring* springA );
	double		getHeight(){ return northPoint->getPos().y - southPoint->getPos().y; }
	double		getWidth(){return 1;}
	//{ return 100.0*(initialWidth-(northPoint->getPos().z-southPoint->getPos().z))/initialWidth;}
	double		getStrain(){ strain = 100*(getHeight()-initialHeight)/initialHeight; return strain;}
	void		setStrain(double newStrain){ strain = newStrain;}
	void		apply( int mode, int nrCorners, Corner* cornerA );
	double		getForceY( int nrCorners, Corner* cornerA);
	void		draw();
private:
	Spring**	pStressSpring;		/**< Spring* array ot pointer to the springs in which the stress is measured. A transsection. */
	int			nrStressSprings;	/**< int number of springs in the transsection. */
	double		stretchVel;			/**< double the velocity [m/s] with which the tissue is stretched/compressed. */
	int			stretching;			/**< int determines in which direction the tissue is stretched: +1: stretch, 0: nothing, -1: compression. */
	double		minStrain;			/**< in % of the original length. */
	double		initialHeight;		/**< */
	double		crossSectionalArea;				/**< the transsectional area of the tissue [m²] */
	double		strain;
	double		northPosY;			// the position of the compressing walls
	double		southPosY;
	Corner*		northPoint;
	Corner*		southPoint;			// to measure the height of the cell in the middle
};

// set the initial values
// ----------------------
void BoundaryConditions::initialize(int nrCorners, Corner* cornerA, Spring* springA)
{
	{ 
		nrStressSprings = 0;
		pStressSpring = new Spring*[nrStressSprings];
	
		// set the north and south corner numbers
		northPoint= &cornerA[nrCorners-1];
		southPoint = &cornerA[0];
		// set up a measure for the width
		initialHeight = this->getHeight();
/*		glColor3f(0,1,0);
		drawDisc(northPoint->getPos(), 0.001);
		glColor3f(0,1,1);
		drawDisc(southPoint->getPos(), 0.001); 
*/
	}

	{ // read from parameterFile.txt
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);

		iF.setScope("-parameters-");
		stretchVel	= iF.getDouble("displacementspeed");
		minStrain	= iF.getDouble("minstrain");
		crossSectionalArea		= iF.getDouble("crosssectionalarea");
		stretching	= -1; //start with compression
	}

	initialHeight = this->getHeight();
	northPosY = northPoint->getPos().y;
	southPosY = southPoint->getPos().y;


	// set the window right
	// the continuation of Container::initialize
	extern double xmax, xmin, ymax, ymin;
	double margin = 2.5*0.5*(-minStrain/100.0)*xmax;
	xmin -= margin;
	xmax += margin;
	ymin -= margin;
	ymax += margin;

}


// apply the BC
// ------------
void BoundaryConditions::apply(int mode, int nrCorners, Corner* cornerA)
{
	if((strain>0        )&&(stretching==+1)){stretching *= -1; /**/}
	if((strain<minStrain)&&(stretching==-1)){stretching *= -1; exit(0);}
	double posY;
	if(mode==1){ // moving
		// apply a sigmoidal acceleration
//		double x = 8.0*((10.0*timeSteps/(double)(1000))-1.0);
//		double boundarySpeed = stretchVel*(1.0/(1.0+(exp(-x))));
		double boundarySpeed = stretchVel;
		northPosY += stretching*boundarySpeed*dt;
		southPosY -= stretching*boundarySpeed*dt;

		
		for(int i=0; i<nrCorners; ++i)
		{
			posY = cornerA[i].getPos().y;
			if(posY>northPosY){
				cornerA[i].setPosY(northPosY);
				cornerA[i].setVel(ZEROTRIPLE);
				cornerA[i].setAcc(ZEROTRIPLE);
				cornerA[i].setBoundary(true);
			}else{
				if(posY<southPosY){
					cornerA[i].setPosY(southPosY);
					cornerA[i].setVel(ZEROTRIPLE);
					cornerA[i].setAcc(ZEROTRIPLE);
					cornerA[i].setBoundary(true);
				}else{
					cornerA[i].setBoundary(false);
				}
			}
		}
	} else { if(mode==0) {// not moving
		for(int i=0; i<nrCorners; ++i)
		{
			posY = cornerA[i].getPos().y;
			if(posY>northPosY){
				cornerA[i].setVel(ZEROTRIPLE);
				cornerA[i].setAcc(ZEROTRIPLE);
				cornerA[i].setBoundary(true);
			}else{
				if(posY<southPosY){
					cornerA[i].setVel(ZEROTRIPLE);
					cornerA[i].setAcc(ZEROTRIPLE);
					cornerA[i].setBoundary(true);
				}else{
					cornerA[i].setBoundary(false);
				}
			}
		}
	} else {cout<<"invalid mode selected for applying BC"<<endl;}}
}

void BoundaryConditions::draw()
{
	glColor4f(0,1,0,0.4);
	double width = 0.0002;
	drawPlane(CTriple(width,0,0), CTriple(0,0,width), CTriple(southPoint->getPos().x-0.5*width,southPosY,southPoint->getPos().z-0.5*width));
	drawPlane(CTriple(width,0,0), CTriple(0,0,width), CTriple(northPoint->getPos().x-0.5*width,northPosY,northPoint->getPos().z-0.5*width));
}

// Calculate the stress in the x-direction in the transsectional area
// ------------------------------------------------------------------
double BoundaryConditions::getForceY(int nrCorners, Corner* cornerA)
{
	// sum up all the forces acting in the y-direction from one hemisphere
	double F = 0;
	int i;
	for(i=0; i<nrCorners/2; ++i)
	{
		F = F + cornerA[i].getF().y;
/*		// for debugging :: draw sampled corners in yellow
		glColor3f(1,1,0);
		drawDisc(cornerA[i].getPos(), 0.00002);
*/	}

	return F; 
}



//////////////////////////////////////////////////////////////////
#endif															//
// end of BoundaryConditions.h									//
//////////////////////////////////////////////////////////////////