//////////////////////////////////////////////////////////////////
#ifndef SPRING_H												//
#define SPRING_H												//
// begin of Spring.h											//
//////////////////////////////////////////////////////////////////

#include "Corner.h"

/**
 * Spring class	
 * models a piece of cell wall, its features are:
 *  - a linear spring connecting two corners 
 *	- restLength keeps track of the length of the spring in rest
 */
class Spring
{	
public:	
	Spring(){}
	~Spring(){}
	void	initialize	( int newNr, Corner* newBegin, Corner* newEnd );
	void	draw		();
	int		getNr		(){ return nr; }
	int		getBeginNr	(){ return pBegin->getNr();}
	int		getEndNr	(){ return pEnd->getNr();}
	int		getBeginNr4SoftTissue(){ return pBegin->getNr4SoftTissue();}
	int		getEndNr4SoftTissue(){ return pEnd->getNr4SoftTissue();}
	void	save		( ofstream outFile );
	void	save4SoftTissue	( ofstream outFile );
	bool	isEqual		( int beginNr, int endNr);
	void	isCB		( bool value){pBegin->isCB(value); pEnd->isCB(value);}
	void	isPB		( bool value){pBegin->isPB(value); pEnd->isPB(value);}
private:
	int		nr;
	Corner*	pBegin;			// position of the begin point
	Corner*	pEnd;			// position of the end point
};

void Spring::initialize(int newNr, Corner* newBegin, Corner* newEnd)
{
	nr		= newNr;
	pBegin	= newBegin;
	pEnd	= newEnd;
}

void Spring::draw()
{
	drawLine(pBegin->getPos(),pEnd->getPos());
}

bool Spring::isEqual(int beginNr, int endNr)
{
	bool	result	= false;
	int		beginN	= pBegin->getNr();
	int		endN	= pEnd->getNr();
	if(((beginNr==beginN)&&(endNr==endN))||((beginNr==endN)&&(endNr==beginN))){
		return true;
	}else{
		return false;
	}
}

void Spring::save(ofstream outFile)
{
	outFile<<nr<<" "<<pBegin->getNr()<<" "<<pEnd->getNr()<<endl; 
}

void Spring::save4SoftTissue(ofstream outFile)
{
	outFile<<nr<<" "<<pBegin->getNr4SoftTissue()<<" "<<pEnd->getNr4SoftTissue()<<endl; 
}





//////////////////////////////////////////////////////////////////
#endif															//
// end of Spring.h												//
//////////////////////////////////////////////////////////////////