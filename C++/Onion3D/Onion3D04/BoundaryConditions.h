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
	void		initialize( Corner* cornerA, Spring* springA, Spring* secSpringA );
	double		getLength(){ return getEastPos() - getWestPos(); }
	double		getWidth(){ return 100.0*(initialWidth-(northPoint->getPos().z-southPoint->getPos().z))/initialWidth;}
	double		getStrain(){ strain = 100*(getLength()-initialLength)/initialLength; return strain;}
	void		setStrain(double newStrain){ strain = newStrain;}
	void		apply( int mode );
	double		getStress();
	void		save( ofstream outFile );
	
private:
	double		getEastPos();
	double		getWestPos();
	Corner**	pEastCorner;		/**< Corner* array of pointers to the corners in the east (right). */
	Corner**	pWestCorner;		/**< Corner* array of pointers to the corners in the west (left). */
	int			nrEastCorners;		/**< int number of corners in the east.*/
	int			nrWestCorners;		/**< int number of corners in the west. */
	Spring**	pStressSpring;		/**< Spring* array ot pointer to the springs in which the stress is measured. A transsection. */
	int			nrStressSprings;	/**< int number of springs in the transsection. */
	Spring**	pStressSecSpring;
	int			nrStressSecSprings;
	double		stretchVel;			/**< double the velocity [m/s] with which the tissue is stretched/compressed. */
	int			stretching;			/**< int determines in which direction the tissue is stretched: +1: stretch, 0: nothing, -1: compression. */
	double		maxStrain;			/**< in % of the original length. */
	double		initialLength;		/**< */
	double		initialWidth;
	double		crossSectionalArea;				/**< the transsectional area of the tissue [m²] */
	double		strain;
	double		eastPosX;
	double		westPosX;
	Corner*		northPoint;
	Corner*		southPoint;			// to measure the width of the tissue in the middle
};

// set the initial values
// ----------------------
void BoundaryConditions::initialize(Corner* cornerA, Spring* springA, Spring* secSpringA)
{
	{ // read from inputFile.txt
		const char* nameIn = "data/inputFile.txt";
		InputFile iF(nameIn);

		iF.setScope("-eastface-");
		nrEastCorners = iF.getInt("nrcorners");
		pEastCorner = new Corner*[nrEastCorners];
		pEastCorner[0] = &cornerA[iF.getInt("numbers")-1];
		int i;
		for(i=1; i<nrEastCorners; ++i)
		{
			pEastCorner[i] = &cornerA[iF.getInt()-1];
		}
		for(i=0; i<nrEastCorners; ++i)
		{
			pEastCorner[i]->setBoundary(true);
		}

		iF.setScope("-westface-");
		nrWestCorners = iF.getInt("nrcorners");
		pWestCorner = new Corner*[nrWestCorners];
		pWestCorner[0] = &cornerA[iF.getInt("numbers")-1];
		for(i=1; i<nrWestCorners; ++i)
		{
			pWestCorner[i] = &cornerA[iF.getInt()-1];
		}
		for(i=0; i<nrWestCorners; ++i)
		{
			pWestCorner[i]->setBoundary(true);
		}

		iF.setScope("-cellboundarycrosssection-");
		nrStressSprings = iF.getInt("nrsprings");
		pStressSpring = new Spring*[nrStressSprings];
		pStressSpring[0] = &springA[iF.getInt("numbers")-1];
		for(i=1; i<nrStressSprings; ++i)
		{
			pStressSpring[i] = &springA[iF.getInt()-1];
		}

		iF.setScope("-planeboundarycrosssection-");
		nrStressSecSprings = iF.getInt("nrsprings");
		pStressSecSpring = new Spring*[nrStressSecSprings];
		pStressSecSpring[0] = &secSpringA[iF.getInt("numbers")-1];
		for(i=1; i<nrStressSecSprings; ++i)
		{
			pStressSecSpring[i] = &secSpringA[iF.getInt()-1];
		}

		// read the north and south corner numbers
		iF.setScope("-northandsouth-");
		northPoint= &cornerA[iF.getInt("northnr")-1];
		southPoint = &cornerA[iF.getInt("southnr")-1];
		// set up a measure for the width
		initialWidth = northPoint->getPos().z-southPoint->getPos().z;
	/*	glColor3f(0,1,0);
		drawDisc(northPoint->getPos(), 0.0001);
		glColor3f(0,1,1);
		drawDisc(southPoint->getPos(), 0.0001); 
	*/
	}

	{ // read from parameterFile.txt
		const char* nameIn = "data/parameterFile.txt";
		InputFile iF(nameIn);

		iF.setScope("-parameters-");
		stretchVel	= iF.getDouble("displacementspeed");
		maxStrain	= iF.getDouble("maxstrain");
		crossSectionalArea		= iF.getDouble("crosssectionalarea");
		stretching	= 1;
	}

	initialLength = this->getLength();
	eastPosX = this->getEastPos();
	westPosX = this->getWestPos();


	// set the window right
	// the continuation of Container::initialize
	extern double xmax, xmin, ymax, ymin;
	double margin = 2.5*0.5*(maxStrain/100.0)*xmax;
	xmin -= margin;
	xmax += margin;
	ymin -= margin;
	ymax += margin;

}

// get the average x-position of the east wall
// -------------------------------------------
double BoundaryConditions::getEastPos()
{
	double eastPosX = 0;
	for(int i=0; i<nrEastCorners; ++i)
	{
		eastPosX += pEastCorner[i]->getPos().x;
	}
	eastPosX /=(double)(nrEastCorners);
	return eastPosX;
}

// get the average x-position of the west wall
// -------------------------------------------
double BoundaryConditions::getWestPos()
{
	double westPosX = 0;
	for(int i=0; i<nrWestCorners; ++i)
	{
		westPosX += pWestCorner[i]->getPos().x;
	}
	westPosX /=(double)(nrWestCorners);
	return westPosX;
}

// apply the BC
// ------------
void BoundaryConditions::apply(int mode)
{
	if((strain<0        )&&(stretching==-1)){stretching *= -1; /*exit(0);*/}
	if((strain>maxStrain)&&(stretching==+1)){stretching *= -1;}
	int i;
	if(mode==1){ // moving
		// apply a sigmoidal acceleration
		double x = 8.0*((10.0*timeSteps/(double)(1000))-1.0);
		stretchVel = 0.5*(1.0/(1.0+(exp(-x))));
		eastPosX += stretching*stretchVel*dt;
		for(i=0; i<nrEastCorners; ++i)
		{
			pEastCorner[i]->setPosX(eastPosX);
			pEastCorner[i]->setVel(ZEROTRIPLE);
			pEastCorner[i]->setAcc(ZEROTRIPLE);
		}
		westPosX -= stretching*stretchVel*dt;
		for(i=0; i<nrWestCorners; ++i)
		{
			pWestCorner[i]->setPosX(westPosX);
			pWestCorner[i]->setVel(ZEROTRIPLE);
			pWestCorner[i]->setAcc(ZEROTRIPLE);
		}
	} else { if(mode==0) {// not moving
		for(i=0; i<nrEastCorners; ++i)
		{
			pEastCorner[i]->setVel(ZEROTRIPLE);
			pEastCorner[i]->setAcc(ZEROTRIPLE);
		}
		for(i=0; i<nrWestCorners; ++i)
		{
			pWestCorner[i]->setVel(ZEROTRIPLE);
			pWestCorner[i]->setAcc(ZEROTRIPLE);
		}
	} else {cout<<"invalid mode selected for applying BC"<<endl;}}
}

// Calculate the stress in the x-direction in the transsectional area
// ------------------------------------------------------------------
double BoundaryConditions::getStress()
{
	double F = 0;
	int i;
	for(i=0; i<nrStressSprings; ++i)
	{
		F = F + pStressSpring[i]->getLoad();
		/* for debugging :: draw stressSprings in red
		glColor3f(1,1,0);
		pStressSpring[i]->draw();*/
	}
	for(i=0; i<nrStressSecSprings; ++i)
	{
		F = F + pStressSecSpring[i]->getLoad();
		/* for debugging :: draw stressSprings in red
		glColor3f(1,1,0);
		pStressSecSpring[i]->draw();	*/
	}
	return F/crossSectionalArea; // divide by the transsectional area in m²
}

// save the boundary conditions
// ----------------------------
void BoundaryConditions::save( ofstream outFile )
{
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrEastCorners<<endl;
	outFile<<"numbers "<<pEastCorner[0]->getNr()<<" ";
	int i;
	for(i=1; i<nrEastCorners; ++i)
	{
		outFile<<pEastCorner[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-eastface-"<<endl;
	outFile<<endl;

	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrWestCorners<<endl;
	outFile<<"numbers "<<pWestCorner[0]->getNr()<<" ";
	for(i=1; i<nrWestCorners; ++i)
	{
		outFile<<pWestCorner[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-westface-"<<endl;
	outFile<<endl;

	outFile<<"-crosssection-"<<endl;
	outFile<<"nrsprings "<<nrStressSprings<<endl;
	outFile<<"numbers "<<pStressSpring[0]->getNr()<<" ";
	for(i=1; i<nrStressSprings; ++i)
	{
		outFile<<pStressSpring[i]->getNr()<<" ";
	}
	outFile<<endl;
	outFile<<"-end-transsection-"<<endl;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of BoundaryConditions.h									//
//////////////////////////////////////////////////////////////////