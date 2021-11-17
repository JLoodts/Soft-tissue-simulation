// The debugger can't handle symbols more than 255 characters long.
// STL often creates symbols longer than that.
// When symbols are longer than 255 characters, the warning is disabled.
// #pragma warning(disable:4786)
 #pragma warning(disable:4715) // not all control paths return a value

/* This program implements an elastic circle composed out of corners which are connected
 * with dampened springs.  An internal pressure is simulated with a simple gamma-gaslaw.
 * 
 * created on 22/01/2003
 * by Jimmy Loodts
 */

#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <fstream.h>	// for file handling
//#include <list>

#include "Main.h"

const	int		nrCllsXL = 3;
const 	int		nrCllsYL = 1;
const	int		nrCllsXT = 24;
const 	int		nrCllsYT = 6;
double tissueDimX = 1440e-6;
double tissueDimY = 120e-6;
const int cellLengthFactor = 4; // cell length is cellLengthFactor times its width (width = dy)

void longitudinal()
{
		const char* nameOut = "long/inputFile.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	const	int		nrCellsX = nrCllsXL;
	const 	int		nrCellsY = nrCllsYL;
	double			dx = tissueDimX/(double)(nrCellsX*cellLengthFactor); 
	double			dy = tissueDimY/(double)(nrCellsY);

	const int		nrDx = nrCellsX*cellLengthFactor + 1;
	const int		nrDy = nrCellsY+1;

	int		nrCorners = nrDx*nrDy;

////////////////////////////////////////////////////////////////////
// get the corners right
////////////////////////////////////////////////////////////////////
	{
		double cornerAx[nrDx][nrDy];
		double cornerAy[nrDx][nrDy];

		double x = 0; double y = 0;
		for(int i=0; i<nrDx; ++i)
		{
			for(int j=0; j<nrDy; ++j)
			{
				cornerAx[i][j] = x;
				cornerAy[i][j] = y;
				y = y + dy;
			}
			x = x + dx;
			y = 0;
		}

		int cornerNr = 1;

		outFile<<"-cornerpoints-"<<endl;
		outFile<<"nrcorners "<<nrCorners<<endl;
		for(i=0; i<nrDx; ++i)
		{
			for(int j=0; j<nrDy; ++j)
			{
				outFile<<cornerNr++<<" "<<cornerAx[i][j]<<" "<<cornerAy[i][j]<<" "<<endl;
			}
		}
		outFile<<"-end-cornerpoints-"<<endl<<endl;
		if(nrCorners!=cornerNr-1)cout<<"Check wether nrcorners is correct!"<<endl;
	}

////////////////////////////////////////////////////////////////////
// get the springs right
////////////////////////////////////////////////////////////////////
	
	int cornerNrA[nrDx][nrDy];	
	
	int cornerNr = 1;
	for(int i=0; i<nrDx; ++i)
	{
		for(int j=0; j<nrDy; ++j)
		{
			cornerNrA[i][j] = cornerNr++;
		}
	}

	const int nrSpringsX = (int)((nrDx-1)*1.5)+1;
	const int nrSpringsY = nrDy;
	int springNrA[nrSpringsX][nrSpringsY];
	for(i=0; i<nrSpringsX; ++i)
	{
		for(int j=0; j<nrSpringsY; ++j)
		{
			springNrA[i][j] = 0;
		}
	}
	int springNr = 1;
	int columnID = 1; // ranges from 1 to cellLengthFactor
	int counter = 0; // counts the number of columns that are written

	for( i=0; i<nrSpringsX; )
	{	
		switch(columnID)
		{
		case 1: // first and last column: all vertical
			{
				for(int j=0; j<nrSpringsY-1; ++j)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 2: // all horizontal
			{
				for(int j=0; j<nrSpringsY; ++j)
				{
					springNrA[i][j] = springNr++;
				}
				counter++;
				if(i<nrSpringsX-2){
					switch(counter)
					{
					case 3: // next one is 3
						{
							columnID = 3;
							i++;
							break;
						}
					case 6: // next one is 4
						{
							columnID = 4;
							i++;
							break;
						}
					default:
						{
							columnID = 2;
							i++;
							break;
						}
					}
				}else{
					columnID = 1; // the last column
					i=nrSpringsX-1;
				}
				break;
			}
		case 3: // even rownumbers vertical
			{
				for(int j=1; j<nrSpringsY-1; j+=2)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 4: // uneven rownumbers vertical
			{
				for(int j=0; j<nrSpringsY-1; j+=2)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter = 1;
				i++;
				break;
			}
		default:
			cout<<"Error in switch(columnID)"<<endl;
			break;
		}

	}
	int nrSprings = springNr-1;

	// lets repeat the same cycles and write our data to file
	outFile<<"-springs-"<<endl;
	outFile<<"nrsprings "<<nrSprings<<endl;

	springNr = 1;
	columnID = 1; // ranges from 1 to cellLengthFactor
	counter = 0; // counts the number of columns that are written

	int ic, jc; ic = jc = 0; // indexes for cornerNrA

	for( i=0; i<nrSpringsX; )
	{	
		switch(columnID)
		{
		case 1: // first and last column: all vertical
			{
				for(int jc=0; jc<nrSpringsY-1; ++jc)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 2: // all horizontal
			{
				for(int jc=0; jc<nrSpringsY; ++jc)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic+1][jc]<<endl;
				}
				counter++;
				if(i<nrSpringsX-2){
					switch(counter)
					{
					case 3: // next one is 3
						{
							columnID = 3;
							i++;
							ic++;
							break;
						}
					case 6: // next one is 4
						{
							columnID = 4;
							i++;
							ic++;
							break;
						}
					default:
						{
							columnID = 2;
							i++;
							ic++;
							break;
						}
					}
				}else{
					columnID = 1; // the last column
					i=nrSpringsX-1;
					ic = nrDx-1;
				}
				break;
			}
		case 3: // even rownumbers vertical
			{
				for(int jc=1; jc<nrSpringsY-1; jc+=2)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter++;
				i++;
				// ic remains the same
				break;
			}
		case 4: // uneven rownumbers vertical
			{
				for(int jc=0; jc<nrSpringsY-1; jc+=2)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter = 1;
				i++;
				// ic remains the same
				break;
			}
		default:
			cout<<"Error in switch(columnID)"<<endl;
			break;
		}

	}
	outFile<<"-end-springs-"<<endl<<endl;

/*
	for(int j=nrSpringsY-1; j>=0; --j)
	{
		for(int i=0; i<nrSpringsX; ++i)
		{
			//cout<<i<<"x"<<j<<" ";
			cout<<springNrA[i][j]<<"  ";
		}
		cout<<endl;
	}
*/	

////////////////////////////////////////////////////////////////////
// get the secundary springs right
////////////////////////////////////////////////////////////////////

	const int nrSecSpringsX = nrDx-1;
	const int nrSecSpringsY = 2*(nrDy-1);
	const int nrSecSprings = nrSecSpringsX*nrSecSpringsY;

	int secSpringNrA[nrSecSpringsX][nrSecSpringsY];
	for(i=0; i<nrSecSpringsX; ++i)
	{
		for(int j=0; j<nrSecSpringsY; ++j)
		{
			secSpringNrA[i][j] = 0;
		}
	}

	outFile<<"-secundarysprings-"<<endl;
	outFile<<"nrsprings "<<nrSecSprings<<endl;

	springNr = 1;

	ic, jc; ic = jc = 0; // indexes for cornerNrA

	for( ic=0; ic<nrSecSpringsX; ++ic)
	{	
		jc = 0;
		for(int j=0; j<nrSecSpringsY; )
		{
			secSpringNrA[ic][j++] = springNr;
			outFile<<springNr++<<" "<<cornerNrA[ic][jc  ]<<" "<<cornerNrA[ic+1][jc+1]<<endl;
			secSpringNrA[ic][j++] = springNr;
			outFile<<springNr++<<" "<<cornerNrA[ic][jc+1]<<" "<<cornerNrA[ic+1][jc  ]<<endl;
			jc++;
		}
	}
	outFile<<"-end-secundarysprings-"<<endl<<endl;

/*	for(j=nrSecSpringsY-1; j>=0; --j)
	{
		for(int i=0; i<nrSecSpringsX; ++i)
		{
			//cout<<i<<"x"<<j<<" ";
			cout<<secSpringNrA[i][j]<<"  ";
		}
		cout<<endl;
	}
*/
////////////////////////////////////////////////////////////////////
// get the cells right
////////////////////////////////////////////////////////////////////
	
	const int nrCells = nrCellsX*nrCellsY+nrCellsY/2;
	int cellNr = 0;
	int icBegin, jcBegin, isBegin, jsBegin, is, js;
	int nrCornersXcell, nrCornersYcell, nrCornersInCell;
	int nrSpringsXcell, nrSpringsYcell, nrSpringsInCell;
	int* cornersOfCell[nrCells];
	int* springsOfCell[nrCells];
	int index;

	// first column of cells
	bool wholeCell = true;
	icBegin = 0; jcBegin = 0; isBegin = 0; jsBegin = 0;
	for(int j=0; j<nrCellsY; ++j)
	{
		if(j%2==0) { // a full cell
			nrCornersXcell = cellLengthFactor+1;
			nrSpringsXcell = cellLengthFactor;
			wholeCell = true;
		} else { // half of a cell
			nrCornersXcell = cellLengthFactor/2+1;
			nrSpringsXcell = cellLengthFactor/2;
			wholeCell = false;
		}
		nrCornersYcell = 2;
		nrSpringsYcell = 1;
		nrCornersInCell = nrCornersXcell*nrCornersYcell;
		nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
		cornersOfCell[cellNr] = new int[nrCornersInCell];
		springsOfCell[cellNr] = new int[nrSpringsInCell];
		cornersOfCell[cellNr][0] = nrCornersInCell;
		springsOfCell[cellNr][0] = nrSpringsInCell;
			
		
		// run a loop such that corners nrs are given in clockwise direction per cell
		index = 1; 
		jc=jcBegin+1;
		for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jc=jcBegin;
		for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jcBegin++;
	
		// run a loop such that spring nrs are given
		index = 1;
		if(wholeCell)
		{
			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[0][js];
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[4][js];
			springsOfCell[cellNr][index++]=springNrA[5][js];
			springsOfCell[cellNr][index++]=springNrA[6][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[4][js];
			springsOfCell[cellNr][index++]=springNrA[5][js];
		} else { // half cell
			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[0][js];
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[3][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
		}
		
		jsBegin++;
		cellNr++;
	}
	isBegin += 3;
	icBegin += 2;

	bool startOnBottom = false;
	// the cells in between
	while(cellNr<nrCellsX*nrCellsY)
	{	 
		if(startOnBottom) { 
			jcBegin = 0; jsBegin = 0;
		} else { 
			jcBegin = 1; jsBegin = 1;
		}
		while(jsBegin<nrCellsY)
		{
			
			nrCornersXcell = cellLengthFactor+1;
			nrSpringsXcell = cellLengthFactor;
			nrCornersYcell = 2;
			nrSpringsYcell = 1;
			nrCornersInCell = nrCornersXcell*nrCornersYcell;
			nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
			cornersOfCell[cellNr] = new int[nrCornersInCell];
			springsOfCell[cellNr] = new int[nrSpringsInCell];
			cornersOfCell[cellNr][0] = nrCornersInCell;
			springsOfCell[cellNr][0] = nrSpringsInCell;
			index = 1;
			
			// run a loop such that corners nrs are given in clockwise direction per cell
			jc=jcBegin+1;
			for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
			{
				cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
			}
			jc=jcBegin;
			for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
			{
				cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
			}
			jcBegin+=2;

				// run a loop such that spring nrs are given
			index = 1;

			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[is+0][js];
			springsOfCell[cellNr][index++]=springNrA[is+1][js];
			springsOfCell[cellNr][index++]=springNrA[is+2][js];
			springsOfCell[cellNr][index++]=springNrA[is+4][js];
			springsOfCell[cellNr][index++]=springNrA[is+5][js];
			springsOfCell[cellNr][index++]=springNrA[is+6][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[is+1][js];
			springsOfCell[cellNr][index++]=springNrA[is+2][js];
			springsOfCell[cellNr][index++]=springNrA[is+4][js];
			springsOfCell[cellNr][index++]=springNrA[is+5][js];
			
			jsBegin+=2;

			cellNr++;
		}
		startOnBottom = !startOnBottom;
		icBegin += 2;
		isBegin += 3;
	}


	// last column of cells
	jcBegin = 1;
	jsBegin = 1;
	for(jsBegin=1; jsBegin<nrCellsY;)
	{
		nrCornersXcell = cellLengthFactor/2+1;
		nrSpringsXcell = cellLengthFactor/2;
		nrCornersYcell = 2;
		nrSpringsYcell = 1;
		nrCornersInCell = nrCornersXcell*nrCornersYcell;
		nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
		cornersOfCell[cellNr] = new int[nrCornersInCell];
		springsOfCell[cellNr] = new int[nrSpringsInCell];
		cornersOfCell[cellNr][0] = nrCornersInCell;
		springsOfCell[cellNr][0] = nrSpringsInCell;
		index = 1;
		
		// run a loop such that corners nrs are given in clockwise direction per cell
		jc=jcBegin+1; js = jsBegin+1;
		for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jc=jcBegin;
		for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jcBegin+=2;

		// run a loop such that spring nrs are given
		index = 1;
		js=jsBegin;
		is=isBegin;
		springsOfCell[cellNr][index++]=springNrA[is][js];
		springsOfCell[cellNr][index++]=springNrA[is+1][js];
		springsOfCell[cellNr][index++]=springNrA[is+2][js];
		springsOfCell[cellNr][index++]=springNrA[is+3][js];
		js++;
		springsOfCell[cellNr][index++]=springNrA[is+1][js];
		springsOfCell[cellNr][index++]=springNrA[is+2][js];
	
		jsBegin+=2;


		cellNr++;
	}

	outFile<<"-configuration-"<<endl;
	outFile<<"nrcells "<<nrCells<<endl;
	for(cellNr=0; cellNr<nrCells; ++cellNr)
	{
		outFile<<cellNr+1<<" ";
		nrCornersInCell = cornersOfCell[cellNr][0];
		nrSpringsInCell = springsOfCell[cellNr][0];
		outFile<<"nrcorners "<<nrCornersInCell<<" numbers ";
		for(int i=1; i<nrCornersInCell+1; ++i){outFile<<cornersOfCell[cellNr][i]<<" ";}
		outFile<<endl;
		outFile<<"   nrsprings "<<nrSpringsInCell<<" numbers ";
		for( i=1; i<nrSpringsInCell+1; ++i){outFile<<springsOfCell[cellNr][i]<<" ";}
		outFile<<endl;
	}
	outFile<<"-end-configuration-"<<endl<<endl;
	
	outFile<<"-northface-"<<endl;
	outFile<<"nrcorners "<<nrDx<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDx; ++i )
	{
		outFile<<cornerNrA[i][nrDy-1]<<" ";
	}
	outFile<<endl<<"-end-northface-"<<endl<<endl;

	
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrDy<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDy; ++i )
	{
		outFile<<cornerNrA[nrDx-1][i]<<" ";
	}
	outFile<<endl<<"-end-eastface-"<<endl<<endl;


	outFile<<"-southface-"<<endl;
	outFile<<"nrcorners "<<nrDx<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDx; ++i )
	{
		outFile<<cornerNrA[i][0]<<" ";
	}
	outFile<<endl<<"-end-southface-"<<endl<<endl;


	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrDy<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDy; ++i )
	{
		outFile<<cornerNrA[0][i]<<" ";
	}
	outFile<<endl<<"-end-westface-"<<endl<<endl;

		
	outFile<<"-crosssection-"<<endl;
	outFile<<"nrsprings "<<nrSpringsY<<endl;
	outFile<<"numbers ";
	int middle = (nrSpringsX/2)-1;
	for(i=0; i<nrSpringsY; ++i)
	{
		outFile<<springNrA[middle][i]<<" ";
	}
	outFile<<endl;
	outFile<<"nrsecundarysprings "<<nrSecSpringsY<<endl;
	outFile<<"secundarynumbers ";
	middle = (nrSecSpringsX/2)-1;
	for(i=0; i<nrSecSpringsY; ++i)
	{
		outFile<<secSpringNrA[middle][i]<<" ";
	}
	outFile<<endl<<"-end-crosssection-"<<endl<<endl;

	outFile.close();
}

///////////////////////////////////////////////////////////
// o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o //
///////////////////////////////////////////////////////////

void transversal()
{
		const char* nameOut = "trans/inputFile.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);
	const	int		nrCellsX = nrCllsYT;
	const 	int		nrCellsY = nrCllsXT;
	double			dx = tissueDimY/(double)(nrCellsX*cellLengthFactor); 
	double			dy = tissueDimX/(double)(nrCellsY);

	const int		nrDx = nrCellsX*cellLengthFactor + 1;
	const int		nrDy = nrCellsY+1;

	int		nrCorners = nrDx*nrDy;

////////////////////////////////////////////////////////////////////
// get the corners right
////////////////////////////////////////////////////////////////////
	{
		double cornerAx[nrDx][nrDy];
		double cornerAy[nrDx][nrDy];

		double x = 0; double y = 0;
		for(int i=0; i<nrDx; ++i)
		{
			for(int j=0; j<nrDy; ++j)
			{
				cornerAx[i][j] = y;
				cornerAy[i][j] = x;
				y = y + dy;
			}
			x = x + dx;
			y = 0;
		}

		int cornerNr = 1;

		outFile<<"-cornerpoints-"<<endl;
		outFile<<"nrcorners "<<nrCorners<<endl;
		for(i=0; i<nrDx; ++i)
		{
			for(int j=0; j<nrDy; ++j)
			{
				outFile<<cornerNr++<<" "<<cornerAx[i][j]<<" "<<cornerAy[i][j]<<" "<<endl;
			}
		}
		outFile<<"-end-cornerpoints-"<<endl<<endl;
		if(nrCorners!=cornerNr-1)cout<<"Check wether nrcorners is correct!"<<endl;
	}

////////////////////////////////////////////////////////////////////
// get the springs right
////////////////////////////////////////////////////////////////////
	
	int cornerNrA[nrDx][nrDy];	
	
	int cornerNr = 1;
	for(int i=0; i<nrDx; ++i)
	{
		for(int j=0; j<nrDy; ++j)
		{
			cornerNrA[i][j] = cornerNr++;
		}
	}

	const int nrSpringsX = (int)((nrDx-1)*1.5)+1;
	const int nrSpringsY = nrDy;
	int springNrA[nrSpringsX][nrSpringsY];
	for(i=0; i<nrSpringsX; ++i)
	{
		for(int j=0; j<nrSpringsY; ++j)
		{
			springNrA[i][j] = 0;
		}
	}
	int springNr = 1;
	int columnID = 1; // ranges from 1 to cellLengthFactor
	int counter = 0; // counts the number of columns that are written

	for( i=0; i<nrSpringsX; )
	{	
		switch(columnID)
		{
		case 1: // first and last column: all vertical
			{
				for(int j=0; j<nrSpringsY-1; ++j)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 2: // all horizontal
			{
				for(int j=0; j<nrSpringsY; ++j)
				{
					springNrA[i][j] = springNr++;
				}
				counter++;
				if(i<nrSpringsX-2){
					switch(counter)
					{
					case 3: // next one is 3
						{
							columnID = 3;
							i++;
							break;
						}
					case 6: // next one is 4
						{
							columnID = 4;
							i++;
							break;
						}
					default:
						{
							columnID = 2;
							i++;
							break;
						}
					}
				}else{
					columnID = 1; // the last column
					i=nrSpringsX-1;
				}
				break;
			}
		case 3: // even rownumbers vertical
			{
				for(int j=1; j<nrSpringsY-1; j+=2)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 4: // uneven rownumbers vertical
			{
				for(int j=0; j<nrSpringsY-1; j+=2)
				{
					springNrA[i][j] = springNr++;
				}
				columnID = 2;
				counter = 1;
				i++;
				break;
			}
		default:
			cout<<"Error in switch(columnID)"<<endl;
			break;
		}

	}
	int nrSprings = springNr-1;

	// lets repeat the same cycles and write our data to file
	outFile<<"-springs-"<<endl;
	outFile<<"nrsprings "<<nrSprings<<endl;

	springNr = 1;
	columnID = 1; // ranges from 1 to cellLengthFactor
	counter = 0; // counts the number of columns that are written

	int ic, jc; ic = jc = 0; // indexes for cornerNrA

	for( i=0; i<nrSpringsX; )
	{	
		switch(columnID)
		{
		case 1: // first and last column: all vertical
			{
				for(int jc=0; jc<nrSpringsY-1; ++jc)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter++;
				i++;
				break;
			}
		case 2: // all horizontal
			{
				for(int jc=0; jc<nrSpringsY; ++jc)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic+1][jc]<<endl;
				}
				counter++;
				if(i<nrSpringsX-2){
					switch(counter)
					{
					case 3: // next one is 3
						{
							columnID = 3;
							i++;
							ic++;
							break;
						}
					case 6: // next one is 4
						{
							columnID = 4;
							i++;
							ic++;
							break;
						}
					default:
						{
							columnID = 2;
							i++;
							ic++;
							break;
						}
					}
				}else{
					columnID = 1; // the last column
					i=nrSpringsX-1;
					ic = nrDx-1;
				}
				break;
			}
		case 3: // even rownumbers vertical
			{
				for(int jc=1; jc<nrSpringsY-1; jc+=2)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter++;
				i++;
				// ic remains the same
				break;
			}
		case 4: // uneven rownumbers vertical
			{
				for(int jc=0; jc<nrSpringsY-1; jc+=2)
				{
					outFile<<springNr++<<" "<<cornerNrA[ic][jc]<<" "<<cornerNrA[ic][jc+1]<<endl;
				}
				columnID = 2;
				counter = 1;
				i++;
				// ic remains the same
				break;
			}
		default:
			cout<<"Error in switch(columnID)"<<endl;
			break;
		}

	}
	outFile<<"-end-springs-"<<endl<<endl;

/*
	for(int j=nrSpringsY-1; j>=0; --j)
	{
		for(int i=0; i<nrSpringsX; ++i)
		{
			//cout<<i<<"x"<<j<<" ";
			cout<<springNrA[i][j]<<"  ";
		}
		cout<<endl;
	}
*/	

////////////////////////////////////////////////////////////////////
// get the secundary springs right
////////////////////////////////////////////////////////////////////

	const int nrSecSpringsX = nrDx-1;
	const int nrSecSpringsY = 2*(nrDy-1);
	const int nrSecSprings = nrSecSpringsX*nrSecSpringsY;

	int secSpringNrA[nrSecSpringsX][nrSecSpringsY];
	for(i=0; i<nrSecSpringsX; ++i)
	{
		for(int j=0; j<nrSecSpringsY; ++j)
		{
			secSpringNrA[i][j] = 0;
		}
	}

	outFile<<"-secundarysprings-"<<endl;
	outFile<<"nrsprings "<<nrSecSprings<<endl;

	springNr = 1;

	ic, jc; ic = jc = 0; // indexes for cornerNrA

	for( ic=0; ic<nrSecSpringsX; ++ic)
	{	
		jc = 0;
		for(int j=0; j<nrSecSpringsY; )
		{
			secSpringNrA[ic][j++] = springNr;
			outFile<<springNr++<<" "<<cornerNrA[ic][jc  ]<<" "<<cornerNrA[ic+1][jc+1]<<endl;
			secSpringNrA[ic][j++] = springNr;
			outFile<<springNr++<<" "<<cornerNrA[ic][jc+1]<<" "<<cornerNrA[ic+1][jc  ]<<endl;
			jc++;
		}
	}
	outFile<<"-end-secundarysprings-"<<endl<<endl;

/*	for(j=nrSecSpringsY-1; j>=0; --j)
	{
		for(int i=0; i<nrSecSpringsX; ++i)
		{
			//cout<<i<<"x"<<j<<" ";
			cout<<secSpringNrA[i][j]<<"  ";
		}
		cout<<endl;
	}
*/
////////////////////////////////////////////////////////////////////
// get the cells right
////////////////////////////////////////////////////////////////////
	
	const int nrCells = nrCellsX*nrCellsY+nrCellsY/2;
	int cellNr = 0;
	int icBegin, jcBegin, isBegin, jsBegin, is, js;
	int nrCornersXcell, nrCornersYcell, nrCornersInCell;
	int nrSpringsXcell, nrSpringsYcell, nrSpringsInCell;
	int* cornersOfCell[nrCells];
	int* springsOfCell[nrCells];
	int index;

	// first column of cells
	bool wholeCell = true;
	icBegin = 0; jcBegin = 0; isBegin = 0; jsBegin = 0;
	for(int j=0; j<nrCellsY; ++j)
	{
		if(j%2==0) { // a full cell
			nrCornersXcell = cellLengthFactor+1;
			nrSpringsXcell = cellLengthFactor;
			wholeCell = true;
		} else { // half of a cell
			nrCornersXcell = cellLengthFactor/2+1;
			nrSpringsXcell = cellLengthFactor/2;
			wholeCell = false;
		}
		nrCornersYcell = 2;
		nrSpringsYcell = 1;
		nrCornersInCell = nrCornersXcell*nrCornersYcell;
		nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
		cornersOfCell[cellNr] = new int[nrCornersInCell];
		springsOfCell[cellNr] = new int[nrSpringsInCell];
		cornersOfCell[cellNr][0] = nrCornersInCell;
		springsOfCell[cellNr][0] = nrSpringsInCell;
			
		
		// run a loop such that corners nrs are given in clockwise direction per cell
		index = 1; 
		jc=jcBegin;
		for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jc=jcBegin+1;
		for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jcBegin++;
	
		// run a loop such that spring nrs are given
		index = 1;
		if(wholeCell)
		{
			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[0][js];
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[4][js];
			springsOfCell[cellNr][index++]=springNrA[5][js];
			springsOfCell[cellNr][index++]=springNrA[6][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[4][js];
			springsOfCell[cellNr][index++]=springNrA[5][js];
		} else { // half cell
			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[0][js];
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
			springsOfCell[cellNr][index++]=springNrA[3][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[1][js];
			springsOfCell[cellNr][index++]=springNrA[2][js];
		}
		
		jsBegin++;
		cellNr++;
	}
	isBegin += 3;
	icBegin += 2;

	bool startOnBottom = false;
	// the cells in between
	while(cellNr<nrCellsX*nrCellsY)
	{	 
		if(startOnBottom) { 
			jcBegin = 0; jsBegin = 0;
		} else { 
			jcBegin = 1; jsBegin = 1;
		}
		while(jsBegin<nrCellsY)
		{
			
			nrCornersXcell = cellLengthFactor+1;
			nrSpringsXcell = cellLengthFactor;
			nrCornersYcell = 2;
			nrSpringsYcell = 1;
			nrCornersInCell = nrCornersXcell*nrCornersYcell;
			nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
			cornersOfCell[cellNr] = new int[nrCornersInCell];
			springsOfCell[cellNr] = new int[nrSpringsInCell];
			cornersOfCell[cellNr][0] = nrCornersInCell;
			springsOfCell[cellNr][0] = nrSpringsInCell;
			index = 1;
			
			// run a loop such that corners nrs are given in clockwise direction per cell
			jc=jcBegin;
			for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
			{
				cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
			}
			jc=jcBegin+1;
			for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
			{
				cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
			}
			jcBegin+=2;

				// run a loop such that spring nrs are given
			index = 1;

			js=jsBegin;
			is=isBegin;
			springsOfCell[cellNr][index++]=springNrA[is+0][js];
			springsOfCell[cellNr][index++]=springNrA[is+1][js];
			springsOfCell[cellNr][index++]=springNrA[is+2][js];
			springsOfCell[cellNr][index++]=springNrA[is+4][js];
			springsOfCell[cellNr][index++]=springNrA[is+5][js];
			springsOfCell[cellNr][index++]=springNrA[is+6][js];
			js++;
			springsOfCell[cellNr][index++]=springNrA[is+1][js];
			springsOfCell[cellNr][index++]=springNrA[is+2][js];
			springsOfCell[cellNr][index++]=springNrA[is+4][js];
			springsOfCell[cellNr][index++]=springNrA[is+5][js];
			
			jsBegin+=2;

			cellNr++;
		}
		startOnBottom = !startOnBottom;
		icBegin += 2;
		isBegin += 3;
	}


	// last column of cells
	jcBegin = 1;
	jsBegin = 1;
	for(jsBegin=1; jsBegin<nrCellsY;)
	{
		nrCornersXcell = cellLengthFactor/2+1;
		nrSpringsXcell = cellLengthFactor/2;
		nrCornersYcell = 2;
		nrSpringsYcell = 1;
		nrCornersInCell = nrCornersXcell*nrCornersYcell;
		nrSpringsInCell = 2*nrSpringsXcell+2*nrSpringsYcell;
		cornersOfCell[cellNr] = new int[nrCornersInCell];
		springsOfCell[cellNr] = new int[nrSpringsInCell];
		cornersOfCell[cellNr][0] = nrCornersInCell;
		springsOfCell[cellNr][0] = nrSpringsInCell;
		index = 1;
		
		// run a loop such that corners nrs are given in clockwise direction per cell
		jc=jcBegin; js = jsBegin+1;
		for(int ic=icBegin; ic<icBegin+nrCornersXcell; ++ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jc=jcBegin+1;
		for(ic=icBegin+nrCornersXcell-1; ic>=icBegin; --ic)
		{
			cornersOfCell[cellNr][index++]=cornerNrA[ic][jc];
		}
		jcBegin+=2;

		// run a loop such that spring nrs are given
		index = 1;
		js=jsBegin;
		is=isBegin;
		springsOfCell[cellNr][index++]=springNrA[is][js];
		springsOfCell[cellNr][index++]=springNrA[is+1][js];
		springsOfCell[cellNr][index++]=springNrA[is+2][js];
		springsOfCell[cellNr][index++]=springNrA[is+3][js];
		js++;
		springsOfCell[cellNr][index++]=springNrA[is+1][js];
		springsOfCell[cellNr][index++]=springNrA[is+2][js];
	
		jsBegin+=2;


		cellNr++;
	}

	outFile<<"-configuration-"<<endl;
	outFile<<"nrcells "<<nrCells<<endl;
	for(cellNr=0; cellNr<nrCells; ++cellNr)
	{
		outFile<<cellNr+1<<" ";
		nrCornersInCell = cornersOfCell[cellNr][0];
		nrSpringsInCell = springsOfCell[cellNr][0];
		outFile<<"nrcorners "<<nrCornersInCell<<" numbers ";
		for(int i=1; i<nrCornersInCell+1; ++i){outFile<<cornersOfCell[cellNr][i]<<" ";}
		outFile<<endl;
		outFile<<"   nrsprings "<<nrSpringsInCell<<" numbers ";
		for( i=1; i<nrSpringsInCell+1; ++i){outFile<<springsOfCell[cellNr][i]<<" ";}
		outFile<<endl;
	}
	outFile<<"-end-configuration-"<<endl<<endl;
	
	outFile<<"-northface-"<<endl;
	outFile<<"nrcorners "<<nrDy<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDy; ++i )
	{
		outFile<<cornerNrA[nrDx-1][i]<<" ";
	}
	outFile<<endl<<"-end-northface-"<<endl<<endl;

	
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrDx<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDx; ++i )
	{
		outFile<<cornerNrA[i][nrDy-1]<<" ";
	}
	outFile<<endl<<"-end-eastface-"<<endl<<endl;


	outFile<<"-southface-"<<endl;
	outFile<<"nrcorners "<<nrDy<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDy; ++i )
	{
		outFile<<cornerNrA[0][i]<<" ";
	}
	outFile<<endl<<"-end-southface-"<<endl<<endl;


	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrDx<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrDx; ++i )
	{
		outFile<<cornerNrA[i][0]<<" ";
	}
	outFile<<endl<<"-end-westface-"<<endl<<endl;

		
	outFile<<"-crosssection-"<<endl;
	int middle = (int)(floor(nrCellsY/2));
	int nr = 0;
	for(i=0; i<nrSpringsX; i++)
	{
		if(springNrA[i][middle]!=0){nr++;}
	}
	outFile<<"nrsprings "<<nr<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrSpringsX; i++)
	{
		if(springNrA[i][middle]!=0)
		{
			outFile<<springNrA[i][middle]<<" ";
		}
	}
	outFile<<endl;
	outFile<<"nrsecundarysprings "<<2*nrSecSpringsX<<endl;
	outFile<<"secundarynumbers ";
	middle = nrSecSpringsY/2;
	for(i=0; i<nrSecSpringsX; ++i)
	{
		outFile<<secSpringNrA[i][middle]<<" ";
		outFile<<secSpringNrA[i][middle+1]<<" ";
	}
	outFile<<endl<<"-end-crosssection-"<<endl<<endl;

	outFile.close();
}

// the main function
// -----------------
int main(int argc, char** argv)
{

		longitudinal();

		transversal();

	return 0;
}