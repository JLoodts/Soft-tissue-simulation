//////////////////////////////////////////////////////////////////
#ifndef CONTAINER_H												//
#define CONTAINER_H												//
// begin of Container.h											//
//////////////////////////////////////////////////////////////////

#include <iostream> 

#include "General.h"
#include "Draw.h"
#include "Corner.h"
#include "Cell.h"
#include "InputFile.h"

extern bool	showPBSprings;
extern bool	showCBSprings;
extern bool	showCells;
extern double DISC_SIZE;
extern GLuint select_buffer[];
extern Pair dimMin, dimMax;  // these values are changed to allow for a small border
extern double tissueThickness;

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
			Container(){}
			Container(Pair newNrGridPoints);
		   ~Container(){};
	void	load				();
	void	draw				();
	void	setSpring			(int beginNr, int endNr);
	void	setSecSpring		(int beginNr, int endNr);
	void	setNorthAndSouth	(int beginNr, int endNr);
	void	setStressSpring		(int beginNr, int endNr);
	void	setSecStressSpring	(int beginNr, int endNr);
	void	setCell				(Pair pos);
	void	setEastCorner		(int beginNr);
	void	setWestCorner		(int beginNr);
	void	removeLastSpring	();
	void	removeLastSecSpring	();
	void	removeLastCell		();
	void	removeLastEastCorner();
	void	removeLastWestCorner();
	void	removeLastStressSpring();
	void	removeLastStressSecSpring();
	void	save				();
	void	save4SoftTissue		();
	int		leftDown			(int x, int y);
	Pair	getCornerPos		(int nr){return cornerA[nr-1].getPos();}
	void	resetNorthAndSouth	(){northNr = southNr = 0;}
private:
	Pair	min;
	Pair	max;
	double	dx;
	double	dy;
	Pair	gridMin;
	Pair	gridD;
	Pair	nrGridPoints;
	int		nrCells;
	int		nrCorners;
	int		nrSprings;
	int		nrSecSprings;
	int		selectedCornerNr;
	int		begin;
	int		end;
	int		nrSelected;
	int		getClosestCorner(Pair pos);
	int		getNextCorner(int currentNr, Pair previousPos, int& springNr);
	int		northNr;
	int		southNr;
	int		nrEastCorners;	
	int		nrWestCorners;	
	int		nrStressSprings;
	int		nrStressSecSprings;
	int*	eastCornerNr;
	int*	westCornerNr;
	int*	stressSpringNr;	
	int*	stressSecSpringNr;
	Cell*	cellA;
	Corner*	cornerA;
	Spring*	springA;
	Spring*	secSpringA;
};



Container::Container(Pair newNrGridPoints)
{
	// Get the corners right
	min				 = dimMin;
	gridMin			 = dimMin;
	max				 = dimMax;
	nrGridPoints	 = newNrGridPoints;
	nrCorners		 = nrGridPoints.x*nrGridPoints.y;
	selectedCornerNr = 0; // 0 means no corner is selected
	begin			 = 0;
	end				 = 0;
	northNr			 = 0;
	southNr			 = 0;
	cornerA			 = new Corner[nrCorners];
	dx				= (max.x - min.x)/(double)(nrGridPoints.x-1);
	dy				= (max.y - min.y)/(double)(nrGridPoints.y-1);
	gridD.x			 = dx;
	gridD.y			 = dy;
	DISC_SIZE		 = dx/10.0;
	Pair	pos		 = min;
	int		counter  = 1;
	int		n		 = 0;
	for(int i=0; i<nrGridPoints.x; ++i){ 
	  for(int j=0; j<nrGridPoints.y; ++j){ 
		  cornerA[n] = Corner(counter,pos);// starts counting from 1 since 0 stands for 'no object picked'
		  counter++; 
		  n++;
		  pos.y += dy;
	  }
	  pos.x += dx;
	  pos.y  = min.y;
	} 

	// Get the cell boundary springs right
	nrSprings	 = 0;
	springA		 = new Spring[2*nrCorners];
	nrSecSprings = 0;
	secSpringA	 = new Spring[2*nrCorners];

	// reserve some space for cells
	nrCells		 = 0;
	cellA		 = new Cell[nrCorners];

	// reserve some space for east and west corners
	nrEastCorners = 0;
	nrWestCorners = 0;
	eastCornerNr  = new int[nrCorners/2];
	westCornerNr  = new int[nrCorners/2];

	// reserve some space for the springs in which the stress will be measured
	nrStressSprings		= 0;
	nrStressSecSprings  = 0;
	stressSpringNr		= new int[nrCorners/2];
	stressSecSpringNr	= new int[nrCorners/2];
	
	// make the boundaries a bit bigger so you see the whole scene
	dimMin.x	 = min.x-0.5*dx;	// these values come from main.cpp such that the main window has a border
	dimMax.x	 = max.x+0.5*dx;
	dimMin.y	 = min.y-0.5*dy;
	dimMax.y	 = max.y+0.5*dy;
	min.x		 = min.x-0.5*dx;	// these are the internal values for the container (for picking)
	max.x		 = max.x+0.5*dx;
	min.y		 = min.y-0.5*dy;
	max.y		 = max.y+0.5*dy;
}

void Container::load()
{
	const char* nameIn = "data/inputFile.txt";
	InputFile iF(nameIn);

	// read the info about the grid
	iF.setScope("-grid-");
	min.x			 = iF.getDouble("minx");
	min.y			 = iF.getDouble("miny");
	gridMin			 = min;
	max.x			 = iF.getDouble("maxx");
	max.y			 = iF.getDouble("maxy");
	nrGridPoints.x	 = iF.getInt("nrgridpointsx");
	nrGridPoints.y	 = iF.getInt("nrgridpointsy");

	// derive the rest
	nrCorners		 = nrGridPoints.x*nrGridPoints.y;
	selectedCornerNr = 0; // 0 means no corner is selected
	begin			 = 0;
	end				 = 0;
	northNr			 = 0;
	southNr			 = 0;
	cornerA			 = new Corner[nrCorners];
	dx				 = (max.x - min.x)/(double)(nrGridPoints.x-1);
	dy				 = (max.y - min.y)/(double)(nrGridPoints.y-1);
	gridD.x			 = dx;
	gridD.y			 = dy;
	DISC_SIZE		 = dx/10.0;
	Pair	pos		 = min;
	int		counter  = 1;
	int		n		 = 0;
	for(int i=0; i<nrGridPoints.x; ++i){ 
	  for(int j=0; j<nrGridPoints.y; ++j){ 
		  cornerA[n] = Corner(counter,pos);// starts counting from 1 since 0 stands for 'no object picked'
		  counter++; 
		  n++;
		  pos.y += dy;
	  }
	  pos.x += dx;
	  pos.y  = min.y;
	} 
	springA		 = new Spring[2*nrCorners];
	secSpringA	 = new Spring[2*nrCorners];
	cellA		 = new Cell[nrCorners];
	eastCornerNr  = new int[nrCorners/2];
	westCornerNr  = new int[nrCorners/2];
	stressSpringNr		= new int[nrCorners/2];
	stressSecSpringNr	= new int[nrCorners/2];


	// read the springs
	iF.setScope("-cellboundarysprings-");
	nrSprings = iF.getInt("nrsprings");
	Corner *begin, *end;
	for(i=0; i<nrSprings; ++i)
	{
		int nr = iF.getInt();
		if(nr!=i+1){cout<<"error in inputFile"<<endl;}
		begin	= &cornerA[iF.getInt()-1];
		end		= &cornerA[iF.getInt()-1];
		springA[i].initialize(nr,begin,end);
		begin->isCB(true);
		end->isCB(true);
	}

	// read the secundary springs
	iF.setScope("-planeboundarysprings-");
	nrSecSprings = iF.getInt("nrsprings");
	for(i=0; i<nrSecSprings; ++i)
	{
		int nr = iF.getInt();
		if(nr!=i+1){cout<<"error in inputFile"<<endl;}
		begin	= &cornerA[iF.getInt()-1];
		end		= &cornerA[iF.getInt()-1];
		secSpringA[i].initialize(nr,begin,end);
		begin->isPB(true);
		end->isPB(true);
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

	// read the east corners
	iF.setScope("-eastface-");
	nrEastCorners = iF.getInt("nrcorners");
	eastCornerNr[0] = iF.getInt("numbers");
	for(i=1; i<nrEastCorners; ++i)
	{
		eastCornerNr[i] = iF.getInt();
	}

	// read the west corners
	iF.setScope("-westface-");
	nrWestCorners = iF.getInt("nrcorners");
	westCornerNr[0] = iF.getInt("numbers");
	for(i=1; i<nrWestCorners; ++i)
	{
		westCornerNr[i] = iF.getInt();
	}

	// read the cell boundary stress springs
	iF.setScope("-cellboundarycrosssection-");
	nrStressSprings = iF.getInt("nrsprings");
	stressSpringNr[0] = iF.getInt("numbers");
	for(i=1; i<nrStressSprings; ++i)
	{
		stressSpringNr[i] = iF.getInt();
	}

	// read the plane boundary stress springs
	iF.setScope("-planeboundarycrosssection-");
	nrStressSecSprings = iF.getInt("nrsprings");
	stressSecSpringNr[0] = iF.getInt("numbers");
	for(i=1; i<nrStressSecSprings; ++i)
	{
		stressSecSpringNr[i] = iF.getInt();
	}

	// read the north and south corner numbers
	iF.setScope("-northandsouth-");
	northNr = iF.getInt("northnr");
	southNr = iF.getInt("southnr");
	
	// make the boundaries a bit bigger so you see the whole scene
	dimMin.x	 = min.x-1.5*dx;	// these values come from main.cpp such that the main window has a border
	dimMax.x	 = max.x+0.5*dx;
	dimMin.y	 = min.y-0.5*dy;
	dimMax.y	 = max.y+0.5*dy;
	min.x		 = min.x-0.5*dx;	// these are the internal values for the container (for picking)
	max.x		 = max.x+0.5*dx;
	min.y		 = min.y-0.5*dy;
	max.y		 = max.y+0.5*dy;
}

// draw the current configuration
// ------------------------------
void Container::draw()
{
	glPushMatrix();
		int i=0;
		if(showCells){
			glColor3f(0.6,0.9,0.6);
			for( i=0; i<nrCells; ++i ){
				cellA[i].drawFull();		
			}
		}
		if(showCBSprings){
			glColor3f(0,0.7,0.3);
			for( i=0; i<nrSprings; ++i ){
				springA[i].draw();		
			}
		}
		if(showPBSprings){
			glColor3f(0.7,0,0.3);
			for( i=0; i<nrSecSprings; ++i ){ 
				secSpringA[i].draw(); 
			}	
		}
		// color is set in corner.draw() according to inUse()
		for( i=0; i<nrCorners; ++i ){
			cornerA[i].draw();		
		}
		if((mode == NORTH_AND_SOUTH)&&(northNr!=0)&&(southNr!=0)){
			glColor3f(1,0.7,0); // orange
			drawDisc(cornerA[southNr-1].getPos(), 1.5*DISC_SIZE); 
			glColor3f(0,0.7,1); // blue
			drawDisc(cornerA[northNr-1].getPos(), 1.5*DISC_SIZE); 
		}
		if(mode == EAST){
			glColor3f(1,0.7,0); // orange
			for(int i=0; i<nrEastCorners; ++i){
				drawDisc(cornerA[eastCornerNr[i]-1].getPos(), 1.5*DISC_SIZE);
			}
		}
		if(mode == WEST){
			glColor3f(1,0.7,0); // orange
			for(int i=0; i<nrWestCorners; ++i){
				drawDisc(cornerA[westCornerNr[i]-1].getPos(), 1.5*DISC_SIZE);
			}
		}
		if(mode == CBTRANS){
			glColor3f(1,0.7,0); // orange
			for(int i=0; i<nrStressSprings; ++i){
				springA[stressSpringNr[i]-1].draw();
			}
		}
		if(mode == PBTRANS){
			glColor3f(1,0.7,0); // orange
			for(int i=0; i<nrStressSecSprings; ++i){
				secSpringA[stressSecSpringNr[i]-1].draw();
			}
		}
	glPopMatrix();
}

int Container::leftDown(int x, int y)
{
	// this is an alternative selection method based on knowledge about the grid
	// transform the pixels into length scales (m)
	double modelMatrix[16];
	double projMatrix[16];
	int    viewport[4];
	double posz; 
	Pair pos;

 	glGetDoublev ( GL_MODELVIEW_MATRIX,  modelMatrix );
	glGetDoublev ( GL_PROJECTION_MATRIX, projMatrix  );
	glGetIntegerv( GL_VIEWPORT,			 viewport	 );           
	gluUnProject ( x,viewport[3]-y,0,modelMatrix,projMatrix,viewport,&pos.x,&pos.y,&posz);
	double	nrPosX	= (pos.x-gridMin.x+0.5*gridD.x)/gridD.x;
	double	nrPosY	= (pos.y-gridMin.y+0.5*gridD.y)/gridD.y;
	int		cellNrX = (int)(int)floor(nrPosX);
	int		cellNrY	= (int)(int)floor(nrPosY);
	int gridPointNr = cellNrY + cellNrX*(int)(nrGridPoints.y) + 1;
	//cout<<"nr of the gridpoint = "<<gridPointNr<<endl;
 
/*	// this is the object selection method from OpenGL
	int gridPointNrBis = 0;
	GLuint    hits;
	GLint     viewport[4];
	GLuint    depth = (GLuint)-1;
	glGetIntegerv	( GL_VIEWPORT, viewport	); // get the current viewport parameters
	glRenderMode	( GL_SELECT				); // set the render mode to selection 
	glInitNames	(						); // clears the name stack
	glPushName	( 0						); // puts 0 on top of the name stack
	//GLuint d[1]; glGetIntegerv(GL_NAME_STACK_DEPTH, d); //-> 3.435.97.3836 names possible
	// setup a picking matrix and render into selection buffer
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
		glLoadIdentity();
		gluPickMatrix(x, viewport[3] - y, 1.0, 1.0, viewport); // must be called BEFORE gluPerspective() or glOrtho()
		gluOrtho2D(min.x, max.x, min.y, max.y);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		this->draw();
		glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	hits = glRenderMode(GL_RENDER);
	if(hits>0) {
		gridPointNrBis = (int)(select_buffer[3]);
	}else{
		gridPointNrBis = 0;
	}
	cout<<gridPointNrBis<<endl;
*/
  return gridPointNr;
  
}

void Container::setSpring(int beginNr, int endNr)
{
	Corner *beginCorner, *endCorner;
	beginCorner	= &cornerA[beginNr-1];
	endCorner	= &cornerA[endNr-1];
	// check wether begin and end are not the same
	bool valid	= true;
	if((beginNr==endNr)||(beginNr==0)||(endNr==0)){
		cout<<"The begin and end points are the same, try again"<<endl;
		valid	= false;
	}
	// check wether this spring isn't already defined
	bool alreadyPresent = false;
	for(int i=0; i<nrSprings; ++i){
		if(springA[i].isEqual(beginNr, endNr)){
			alreadyPresent = true;
			cout<<"this spring is already present"<<endl;
		}
	}
	if((valid)&&(!alreadyPresent)) {
		springA[nrSprings].initialize(nrSprings+1,beginCorner,endCorner);
		springA[nrSprings].isCB(true);
		nrSprings++;
	}
}

void Container::setSecSpring(int beginNr, int endNr)
{
	Corner *beginCorner, *endCorner;
	beginCorner	= &cornerA[beginNr-1];
	endCorner	= &cornerA[endNr-1];
	// check wether begin and end are not the same
	bool valid	= true;
	if((beginNr==endNr)||(beginNr==0)||(endNr==0)){
		cout<<"The begin and end points are the same, try again"<<endl;
		valid	= false;
	}
	// check wether this spring isn't already defined
	bool alreadyPresent = false;
	for(int i=0; i<nrSecSprings; ++i){
		if(secSpringA[i].isEqual(beginNr, endNr)){
			alreadyPresent = true;
			cout<<"this spring is already present"<<endl;
		}
	}
	if((valid)&&(!alreadyPresent)) {
		secSpringA[nrSecSprings].initialize(nrSecSprings+1,beginCorner,endCorner);
		secSpringA[nrSecSprings].isPB(true);
		nrSecSprings++;
		cout<<"sec spring created between "<<beginNr<<" and "<<endNr<<endl;
	}
}

void Container::setNorthAndSouth(int beginNr, int endNr)
{
	bool valid	= true;
	if((beginNr==endNr)||(beginNr==0)||(endNr==0)){
		cout<<"The begin and end points are the same, try again"<<endl;
		valid	= false;
	}
	if((!cornerA[beginNr-1].isCB())||(!cornerA[endNr-1].isCB())){
		cout<<"One of these points is not a cell boundary point, try again"<<endl;
		valid = false;
	}
	if(valid){
		southNr = beginNr;
		northNr = endNr;
	}
}

void Container::setEastCorner(int beginNr)
{
	bool valid = true;
	if(beginNr==0){
		valid = false;
	}
	if(!cornerA[beginNr-1].isCB()){
		cout<<"This point is not a cell boundary point, try again"<<endl;
		valid = false;
	}
	if(valid){
		eastCornerNr[nrEastCorners] = beginNr;
		nrEastCorners++;
	}
}

void Container::setWestCorner(int beginNr)
{
	bool valid = true;
	if(beginNr==0){
		valid = false;
	}
	if(!cornerA[beginNr-1].isCB()){
		cout<<"This point is not a cell boundary point, try again"<<endl;
		valid = false;
	}
	if(valid){
		westCornerNr[nrWestCorners] = beginNr;
		nrWestCorners++;
	}
}

void Container::setStressSpring(int beginNr, int endNr)
{
	bool valid	= true;
	if((beginNr==endNr)||(beginNr==0)||(endNr==0)){
		cout<<"The begin and end points are the same, try again"<<endl;
		valid	= false;
	}
	int nrSelectedSpring = 0;
	for(int i=0; i<nrSprings; ++i){
		if(springA[i].isEqual(beginNr,endNr)){
			nrSelectedSpring = i+1;
		}
	}
	if((valid)&&(nrSelectedSpring)){
		stressSpringNr[nrStressSprings] = nrSelectedSpring;
		nrStressSprings++;
	}
}

void Container::setSecStressSpring(int beginNr, int endNr)
{
	bool valid	= true;
	if((beginNr==endNr)||(beginNr==0)||(endNr==0)){
		cout<<"The begin and end points are the same, try again"<<endl;
		valid	= false;
	}
	int nrSelectedSecSpring = 0;
	for(int i=0; i<nrSecSprings; ++i){
		if(secSpringA[i].isEqual(beginNr,endNr)){
			nrSelectedSecSpring = i+1;
		}
	}
	if((valid)&&(nrSelectedSecSpring)){
		stressSecSpringNr[nrStressSecSprings] = nrSelectedSecSpring;
		nrStressSecSprings++;
	}
}

void Container::removeLastSpring()
{
	if(nrSprings>0){
		springA[nrSprings-1].isCB(false);				
		nrSprings--;
		for(int i=0; i<nrSprings; ++i ){
			springA[i].isCB(true);		
		}
	}
}

void Container::removeLastSecSpring()
{
	if(nrSecSprings>0){
		secSpringA[nrSecSprings-1].isPB(false);				
		nrSecSprings--;
		for(int i=0; i<nrSecSprings; ++i ){
			secSpringA[i].isPB(true);		
		}
	}
}

void Container::removeLastCell()
{
	if(nrCells>0){
		nrCells--;
	}
}

void Container::removeLastEastCorner()
{
	if(nrEastCorners>0){
		nrEastCorners--;
	}
}

void Container::removeLastWestCorner()
{
	if(nrWestCorners>0){
		nrWestCorners--;
	}
}

void Container::removeLastStressSpring()
{
	if(nrStressSprings>0){
		nrStressSprings--;
	}
}

void Container::removeLastStressSecSpring()
{
	if(nrStressSecSprings>0){
		nrStressSecSprings--;
	}
}

int	Container::getClosestCorner(Pair pos)
{
	double	nrPosX	= (pos.x-gridMin.x+0.5*gridD.x)/gridD.x;
	double	nrPosY	= (pos.y-gridMin.y+0.5*gridD.y)/gridD.y;
	int		cellNrX = (int)(int)floor(nrPosX);
	int		cellNrY	= (int)(int)floor(nrPosY);
	int closestNr = cellNrY + cellNrX*(int)(nrGridPoints.y) + 1;

	/*// can be made much more efficient by including knowledge about the grid!
	int closestNr = 0;
	double distance = (max.x-min.x) + (max.y-min.y);
	double proximity = distance;
	for(int i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()){
			distance = Distance(pos,cornerA[i].getPos());
			if(distance<proximity){
				proximity = distance;
				closestNr = cornerA[i].getNr();
			}
		}
	}
	*/
	return closestNr;
}

int Container::getNextCorner(int currentNr, Pair previousPos, int& springNr)
{
	double testAngle, angle = 360; // degrees
	springNr = 0;
	int testNr, nextNr;
	Pair currentPos = cornerA[currentNr-1].getPos();
	for(int i=0; i<nrSprings; ++i)
	{
		if(springA[i].getBeginNr()==currentNr){
			testNr = springA[i].getEndNr();
			testAngle = AngleDeg(previousPos-currentPos,cornerA[testNr-1].getPos()-currentPos);
			if(testAngle<=0){testAngle+=360;}
			if(testAngle<angle){
				nextNr = testNr;
				angle = testAngle;
				springNr = i+1;
			}
		}else{
			if(springA[i].getEndNr()==currentNr){
				testNr = springA[i].getBeginNr();
				testAngle = AngleDeg(previousPos-currentPos,cornerA[testNr-1].getPos()-currentPos);
				if(testAngle<=0){testAngle+=360;}
				if(testAngle<angle){
					nextNr = testNr;
					angle = testAngle;
					springNr = i+1;
				}
			}
		}
	}
	return nextNr;
}

void Container::setCell(Pair pos)
{
	int nextNr, currentNr, initialNr, springNr, newNrSprings = 0, newNrCorners = 0;
	Pair previousPos;
	// first run is just to get the number of springs and corners right...
	initialNr		= currentNr = getClosestCorner(pos);
	previousPos		= pos;
	nextNr			= getNextCorner(currentNr,previousPos,springNr);
	while((nextNr!=initialNr)&&(newNrCorners<nrCorners))
	{
		nextNr = getNextCorner(currentNr,previousPos,springNr);
		previousPos = cornerA[currentNr-1].getPos(); 
		currentNr = nextNr;
		newNrCorners++;
	}
	newNrSprings	= newNrCorners;
	int* cornerNrs;
	int* springNrs;
	cornerNrs		= new int[newNrCorners];
	springNrs		= new int[newNrSprings];
	// this is the real run...
	initialNr		= currentNr = getClosestCorner(pos);
	previousPos		= pos;
	nextNr			= getNextCorner(currentNr,previousPos,springNr);
	int nrC			= 0;
//	cout<<"--------selected------"<<endl;
	while((nextNr!=initialNr)&&(nrC<=newNrCorners))
	{	
		nextNr = getNextCorner(currentNr,previousPos,springNr);
		previousPos = cornerA[currentNr-1].getPos(); 
		currentNr = nextNr;
		cornerNrs[nrC] = currentNr;
		springNrs[nrC] = springNr;	
//		cout<<"corner "<<currentNr<<"   spring "<<springNr<<endl;
		nrC++;
	}
	if(newNrCorners>=3){
//		cout<<"--------approved------ NR "<<nrCells<<endl;
		cellA[nrCells].initialize(nrCells+1, newNrCorners, cornerNrs, cornerA, 
									         newNrSprings, springNrs, springA );
		nrCells++;
	}
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
	
	int i, count = 0;
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			count++;
		}
	}
	outFile<<"nrcorners "<<count<<endl;
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			cornerA[i].save(outFile);
		}
	}
	outFile<<"-end-cornerpoints-"<<endl;
	outFile<<endl;
	
	// save the cell boundary springs
	outFile<<"-cellboundarysprings-"<<endl;
	outFile<<"nrsprings "<<nrSprings<<endl;
	for(i=0; i<nrSprings; ++i)
	{
		springA[i].save(outFile);
	}
	outFile<<"-end-cellboundarysprings-"<<endl;
	outFile<<endl;

	// save the plane boundary springs
	outFile<<"-planeboundarysprings-"<<endl;
	outFile<<"nrsprings "<<nrSecSprings<<endl;
	for(i=0; i<nrSecSprings; ++i)
	{
		secSpringA[i].save(outFile);
	}
	outFile<<"-end-planeboundarysprings-"<<endl;
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

	// save the boundary conditions
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrEastCorners<<endl;
	outFile<<"numbers "<<eastCornerNr[0]<<" ";
	for(i=1; i<nrEastCorners; ++i)
	{
		outFile<<eastCornerNr[i]<<" ";
	}
	outFile<<endl;
	outFile<<"-end-eastface-"<<endl;
	outFile<<endl;

	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrWestCorners<<endl;
	outFile<<"numbers "<<westCornerNr[0]<<" ";
	for(i=1; i<nrWestCorners; ++i)
	{
		outFile<<westCornerNr[i]<<" ";
	}
	outFile<<endl;
	outFile<<"-end-westface-"<<endl;
	outFile<<endl;

	outFile<<"-cellboundarycrosssection-"<<endl;
	outFile<<"nrsprings "<<nrStressSprings<<endl;
	outFile<<"numbers "<<stressSpringNr[0]<<" ";
	for(i=1; i<nrStressSprings; ++i)
	{
		outFile<<stressSpringNr[i]<<" ";
	}
	outFile<<endl;
	outFile<<"-end-cellboundarycrosssection-"<<endl;
	outFile<<endl;
	
	outFile<<"-planeboundarycrosssection-"<<endl;
	outFile<<"nrsprings "<<nrStressSecSprings<<endl;
	outFile<<"numbers "<<stressSecSpringNr[0]<<" ";
	for(i=1; i<nrStressSecSprings; ++i)
	{
		outFile<<stressSecSpringNr[i]<<" ";
	}
	outFile<<endl;
	outFile<<"-end-planeboundarycrosssection-"<<endl;
	outFile<<endl;

	outFile<<"-northandsouth-"<<endl;
	outFile<<"northnr "<<northNr<<endl;
	outFile<<"southnr "<<southNr<<endl;
	outFile<<"-end-northandsouth-"<<endl;
	outFile<<endl;

	outFile<<"-grid-"<<endl;
	outFile<<"minx "<<min.x+0.5*dx<<endl;
	outFile<<"miny "<<min.y+0.5*dy<<endl;
	outFile<<"maxx "<<max.x-0.5*dx<<endl;
	outFile<<"maxy "<<max.y-0.5*dy<<endl;
	outFile<<"nrgridpointsx "<<nrGridPoints.x<<endl;
	outFile<<"nrgridpointsy "<<nrGridPoints.y<<endl;
	outFile<<"-end-grid-"<<endl;

}

// save everything so it can be read into SoftTissue
// -------------------------------------------------
void Container::save4SoftTissue()
{
	const char* nameOut = "data/inputFile.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	if (! outFile){}

	// save the cornerpoints
	outFile<<"-cornerpoints-"<<endl;
	
	// renumber the corners by setting their nr4SoftTissue
	int i, nrActiveCorners = 0;
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			cornerA[i].setNr4SoftTissue(nrActiveCorners+1);
			nrActiveCorners++;
		}
	}

	// every corner is saved twice once with z = 0 and a second time with z = d
	outFile<<"nrcorners "<<2*nrActiveCorners<<endl;
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			cornerA[i].save4SoftTissue(outFile,0);
		}
	}
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			outFile<<cornerA[i].getNr4SoftTissue()+nrActiveCorners<<" "<<cornerA[i].getPos().x<<" "<<tissueThickness<<" "<<cornerA[i].getPos().y<<endl;
		}
	}
	outFile<<"-end-cornerpoints-"<<endl;
	outFile<<endl;
	
	// save the cell boundary springs
	outFile<<"-cellboundarysprings-"<<endl;
	outFile<<"nrsprings "<<4*nrSprings+nrActiveCorners<<endl;
	for(i=0; i<nrSprings; ++i)
	{
		springA[i].save4SoftTissue(outFile);
	}
	for(i=0; i<nrSprings; ++i)
	{
		outFile<<springA[i].getNr()+nrSprings<<" "<<springA[i].getBeginNr4SoftTissue()+nrActiveCorners<<" "<<springA[i].getEndNr4SoftTissue()+nrActiveCorners<<endl; 
	}
	int counter = 0;
	// these are the springs parallel to the y-axis
	for(i=0; i<nrCorners; ++i)
	{
		if(cornerA[i].isCB()||cornerA[i].isPB()){
			counter++;
			outFile<<2*nrSprings+counter<<" "<<cornerA[i].getNr4SoftTissue()<<" "<<cornerA[i].getNr4SoftTissue()+nrActiveCorners<<endl;
		}
	}
	// these springs are the diagonals on the side-faces
	for(i=0; i<nrSprings; ++i)
	{
		counter++;
		outFile<<2*nrSprings+counter<<" "<<springA[i].getBeginNr4SoftTissue()<<" "<<springA[i].getEndNr4SoftTissue()+nrActiveCorners<<endl; 
		counter++;
		outFile<<2*nrSprings+counter<<" "<<springA[i].getBeginNr4SoftTissue()+nrActiveCorners<<" "<<springA[i].getEndNr4SoftTissue()<<endl; 
	}
	outFile<<"-end-cellboundarysprings-"<<endl;
	outFile<<endl;

	// not needed in 3D
	// save the plane boundary springs
	outFile<<"-planeboundarysprings-"<<endl;
	outFile<<"nrsprings "<<nrSecSprings<<endl;
	for(i=0; i<nrSecSprings; ++i)
	{
		secSpringA[i].save4SoftTissue(outFile);
	}
	outFile<<"-end-planeboundarysprings-"<<endl;
	outFile<<endl;

	// save the configuration
	outFile<<"-configuration-"<<endl;
	outFile<<"nrcells "<<nrCells<<endl;
	for(i=0; i<nrCells; ++i)
	{
		cellA[i].save4SoftTissue(outFile,nrActiveCorners);
	}
	outFile<<"-end-configuration-"<<endl;
	outFile<<endl;

	// save the boundary conditions
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<2*nrEastCorners<<endl;
	outFile<<"numbers "<<cornerA[eastCornerNr[0]-1].getNr4SoftTissue()<<" ";
	for(i=1; i<nrEastCorners; ++i)
	{
		outFile<<cornerA[eastCornerNr[i]-1].getNr4SoftTissue()<<" ";
	}
	for(i=0; i<nrEastCorners; ++i)
	{
		outFile<<cornerA[eastCornerNr[i]-1].getNr4SoftTissue()+nrActiveCorners<<" ";
	}
	outFile<<endl;
	outFile<<"-end-eastface-"<<endl;
	outFile<<endl;

	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<2*nrWestCorners<<endl;
	outFile<<"numbers "<<cornerA[westCornerNr[0]-1].getNr4SoftTissue()<<" ";
	for(i=1; i<nrWestCorners; ++i)
	{
		outFile<<cornerA[westCornerNr[i]-1].getNr4SoftTissue()<<" ";
	}
	for(i=0; i<nrWestCorners; ++i)
	{
		outFile<<cornerA[westCornerNr[i]-1].getNr4SoftTissue()+nrActiveCorners<<" ";
	}
	outFile<<endl;
	outFile<<"-end-westface-"<<endl;
	outFile<<endl;

	outFile<<"-cellboundarycrosssection-"<<endl;
	outFile<<"nrsprings "<<2*nrStressSprings<<endl;
	outFile<<"numbers "<<stressSpringNr[0]<<" ";
	for(i=1; i<nrStressSprings; ++i)
	{
		outFile<<stressSpringNr[i]<<" ";
	}
	for(i=0; i<nrStressSprings; ++i)
	{
		outFile<<stressSpringNr[i]+nrSprings<<" ";
	}
	outFile<<endl;
	outFile<<"-end-cellboundarycrosssection-"<<endl;
	outFile<<endl;
	
	// not needed in 3D
	outFile<<"-planeboundarycrosssection-"<<endl;
	outFile<<"nrsprings "<<nrStressSecSprings<<endl;
	outFile<<"numbers "<<stressSecSpringNr[0]<<" ";
	for(i=1; i<nrStressSecSprings; ++i)
	{
		outFile<<stressSecSpringNr[i]<<" ";
	}
	outFile<<endl;
	outFile<<"-end-planeboundarycrosssection-"<<endl;
	outFile<<endl;

	outFile<<"-northandsouth-"<<endl;
	outFile<<"northnr "<<cornerA[northNr-1].getNr4SoftTissue()<<endl;
	outFile<<"southnr "<<cornerA[southNr-1].getNr4SoftTissue()<<endl;
	outFile<<"-end-northandsouth-"<<endl;
	outFile<<endl;

	outFile<<"-grid-"<<endl;
	outFile<<"minx "<<min.x+0.5*dx<<endl;
	outFile<<"miny "<<min.y+0.5*dy<<endl;
	outFile<<"maxx "<<max.x-0.5*dx<<endl;
	outFile<<"maxy "<<max.y-0.5*dy<<endl;
	outFile<<"nrgridpointsx "<<nrGridPoints.x<<endl;
	outFile<<"nrgridpointsy "<<nrGridPoints.y<<endl;
	outFile<<"-end-grid-"<<endl;

}
//////////////////////////////////////////////////////////////////
#endif															//
// end of Container.h											//
//////////////////////////////////////////////////////////////////