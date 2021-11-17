//////////////////////////////////////////////////////////////////
#ifndef CORNER_H												//
#define CORNER_H												//
// begin of Corner.h											//
//////////////////////////////////////////////////////////////////

#include "General.h"
#include <fstream.h>	// for file handling
#include "Container.h"

extern double	DISC_SIZE;

/**
 *  The Corner class. 
 *	
 */
class Corner
{	
public:
	Corner(){}
	Corner(int newNr, Pair newPos);
	~Corner(){}
	int		getNr	(					){return nr;}
	Pair	getPos	(					){return pos;}
	void	draw	();
	void	isCB	( bool value		){thisIsCB = value;}
	void	isPB	( bool value		){thisIsPB = value;}
	bool	isCB	(					){return thisIsCB;}
	bool	isPB	(					){return thisIsPB;}
	void	save	( ofstream outFile  ){ outFile<<nr<<" "<<pos.x<<" "<<pos.y<<endl; }
	void	save4SoftTissue	( ofstream outFile  );
	void	setNr4SoftTissue( int value){nr4SoftTissue = value;}
	int		getNr4SoftTissue(){return nr4SoftTissue;}
private:
	int		nr;				/**< from 1 to nrGridPoints. */
	int		nr4SoftTissue;
	Pair	pos;			/**< Pair Position [m] of the corner. */
	bool	thisIsCB;
	bool	thisIsPB;
};

Corner::Corner(int newNr, Pair newPos)
{
	nr = newNr; 
	nr4SoftTissue = 0; 
	pos = newPos; 
	thisIsCB = false; 
	thisIsPB = false;
}

void Corner::draw()
{
	glPushMatrix();
		glLoadName(nr);
		if(thisIsCB){
			glColor3f(0.6,0.6,1);
		}else{
			if(thisIsPB){
				glColor3f(0,1,0);
			}else{
			glColor3f(0.6,0.6,0.6);
			}
		}
		drawDisc(pos,DISC_SIZE);
	glPopMatrix();
}

void Corner::save4SoftTissue( ofstream outFile  )
{
	outFile<<nr4SoftTissue<<" "<<pos.x<<" "<<pos.y<<endl;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of Corner.h												//
//////////////////////////////////////////////////////////////////