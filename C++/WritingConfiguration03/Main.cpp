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

const	int		nrCllsXL = 37;
const 	int		nrCllsYL = 25;
const	int		nrCllsXT = 7;
const 	int		nrCllsYT = 151;
double cllDimX = 600e-6;
double cllDimY = 120e-6;

void longitudinal()
{
		const char* nameOut = "data/inputFileL.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);

	double cellDimX = cllDimX;
	double cellDimY = cllDimY;
	int dxWidth = 6;
	double dx = cellDimX/(double)(dxWidth+2); 
	double dy = cellDimY/(double)(2);
const	int		nrCellsX = nrCllsXL;
const 	int		nrCellsY = nrCllsYL;
	int		nrDX = (dxWidth+1)*nrCellsX+1;
	int		nrDY = 2*nrCellsY;

	int		bigColumnDY = nrDY+1;
	int		smallColumnDY = nrDY;
	int		columnCode = 1;

	int cornerNr = 0;
	double x,y;
	outFile<<"-cornerpoints-"<<endl;
	for(int ix=0; ix<=nrDX; )
	{
		x = ix*dx;
		if((columnCode==1)||(columnCode==4))
		{
			for(int iy=1; iy<=smallColumnDY; iy += 2)
			{
				 y = iy*dy;
				 outFile<<++cornerNr<<" "<<x<<" "<<y<<endl;
			}
			if(columnCode==1){columnCode=2; ix+=1;}
			else{/*(columnCode==4)*/columnCode=1; ix+=dxWidth;}
		} else { // (columnCode==2)||(columnCode==3)
			for(int iy=0; iy<=bigColumnDY; iy += 2)
			{
				 y = iy*dy;
				 outFile<<++cornerNr<<" "<<x<<" "<<y<<endl;
			}
			if(columnCode==2){columnCode=3; ix+=dxWidth;}
			else{/*(columnCode==3)*/columnCode=4; ix+=1;}
		}
	}
	int nrCorners = cornerNr;
	outFile<<"nrcorners "<<nrCorners<<endl;
	outFile<<"-end-cornerpoints-"<<endl<<endl;

////////////////////////////////////////////////////////////////////
// get the springs right
////////////////////////////////////////////////////////////////////
	int springNr = 0;
	columnCode = 1;
	int rightBegin = 1;
	int left;
	int	right;
	outFile<<"-springs-"<<endl;
	int i;
	for(ix=0; ix<nrDX;  )
	{
		if(columnCode==1)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY;
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<left  <<" "<<right++<<endl;
				 outFile<<++springNr<<" "<<left++<<" "<<right  <<endl;
			}
			{columnCode=2; ix+=1;}
		}else{
		if(columnCode==2)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY+1;
			right = rightBegin;
			for(i=0; i<nrCellsY+1; ++i)
			{
				 outFile<<++springNr<<" "<<left++<<" "<<right++<<endl;
			}
			{columnCode=3; ix+=dxWidth;}
		}else{
		if(columnCode==3)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY+1;
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<left++<<" "<<right  <<endl;
				 outFile<<++springNr<<" "<<left  <<" "<<right++<<endl;
			}
			{columnCode=4; ix+=1;}
		}else{
		if(columnCode==4)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY; // one less than the others!
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<left++<<" "<<right++<<endl;
			}
			{columnCode=1; ix+=dxWidth;}
		}
		} // end if(columnCode==1)else
		} // end if(columnCode==2)else
		} // end if(columnCode==3)else
	}
	int nrSprings = springNr;
	outFile<<"nrsprings "<<nrSprings<<endl;
	outFile<<"-end-springs-"<<endl<<endl;

////////////////////////////////////////////////////////////////////
// get the secundary springs right
////////////////////////////////////////////////////////////////////
	springNr = 0;
	int wi, w, sw, nw, se, ne, e;
	bool thickColumn = true;
	wi = w = 1;
	sw = w + nrCellsY;
	nw = sw + 1;
	se = sw + nrCellsY + 1;
	ne = se + 1;
	e = se + nrCellsY + 1;

	outFile<<"-secundarysprings-"<<endl;
	int cellColumnX;
	int cellRowY = 0;
	for(cellColumnX=0; cellColumnX<nrCellsX;  )
	{
		while(cellRowY<nrCellsY)
		{
//			outFile<<++springNr<<" "<<sw<<" "<<e <<endl;
//			outFile<<++springNr<<" "<<se<<" "<<w <<endl;
			outFile<<++springNr<<" "<<sw<<" "<<nw<<endl;
			outFile<<++springNr<<" "<<se<<" "<<ne<<endl;
//			outFile<<++springNr<<" "<<w <<" "<<ne<<endl;
//			outFile<<++springNr<<" "<<e <<" "<<nw<<endl;

			sw = nw;
			se = ne;
			w++;
			e++;
			nw++;
			ne++;
			cellRowY++;
			
		}
		

		if(thickColumn){
			wi += 2*nrCellsY + 2;
			w = wi;
			sw = w + nrCellsY;
			nw = sw + 1;
			se = sw + nrCellsY;
			ne = se + 1;
			e  = se + nrCellsY + 1;
			cellRowY = 1;
		} else {
			wi += 2*nrCellsY;
			w = wi;
			sw = w + nrCellsY;
			nw = sw + 1;
			se = sw + nrCellsY + 1;
			ne = se + 1;
			e  = se + nrCellsY + 1;
			cellRowY = 0;
		}
		cellColumnX++;
		thickColumn = !thickColumn;
	}
	int nrSecundarySprings = springNr;
	outFile<<"nrsprings "<<nrSecundarySprings<<endl;
	outFile<<"-end-secundarysprings-"<<endl<<endl;
	
////////////////////////////////////////////////////////////////////
// get the cells right
////////////////////////////////////////////////////////////////////

	int wallNr = 1;
	int wallColumn=1;
	int wall[2*nrCellsY][2*nrCellsX+1];
	for(int j=0; j<2*nrCellsX+1; ++j){
		if((wallColumn==1)||(wallColumn==3))
		{
			for(i=0; i<2*nrCellsY; ++i)
			{
				wall[i][j] = wallNr++;
			}
			wallColumn++;
		}else
		{
			if(wallColumn==2)
			{
				for(i=0; i<2*nrCellsY; ++i)
				{
					if(i%2!=0)
					{
						wall[i][j] = wallNr;
					}else
					{
						wall[i][j] = wallNr++;
					}	
				}
				wallNr++;
				wallColumn = 3;
			}else
			{
				if(wallColumn==4)
				{
					wall[0][j] = 0;
					for(i=1; i<2*nrCellsY-1; ++i)
					{
						if(i%2==0)
						{
							wall[i][j] = wallNr;
						}else
						{
							wall[i][j] = wallNr++;
						}	
					}
					wallNr++;
					wall[2*nrCellsY-1][j] = 0;
					wallColumn = 1;
				}
			}
		}
	}

	outFile<<"-configuration-"<<endl;
	// first row of cells
	int cellNr=0;
	cornerNr = 1;

	for(i=0; i<nrCellsX; ++i)
	{
		if(i%2==0) // a big column of nrCellsY
		{
			for(int j=0; j<nrCellsY; ++j)
			{
				outFile<<++cellNr
					<<" nrcorners "<<6
					<<" numbers "
					<<cornerNr<<" "<<cornerNr+nrCellsY+1<<" "<<cornerNr+2*nrCellsY+2<<" "<<cornerNr+3*nrCellsY+2<<" "<<cornerNr+2*nrCellsY+1<<" "<<cornerNr+nrCellsY<<endl
					<<"  nrsprings "<<6
					<<" numbers "
					<<wall[2*j][2*i]<<" "<<wall[2*j+1][2*i]<<" "<<wall[2*j+1][2*i+1]<<" "<<wall[2*j+1][2*i+2]<<" "<<wall[2*j][2*i+2]<<" "<<wall[2*j][2*i+1]<<endl;
				cornerNr++;
			}
			cornerNr += nrCellsY+2;
			wallNr += nrCellsY+2;
		}else{ // a small column of nrCellsY-1
			for(int j=0; j<nrCellsY-1; ++j)
			{
				outFile<<++cellNr
					<<" nrcorners "<<6
					<<" numbers "
					<<cornerNr<<" "<<cornerNr+nrCellsY+1<<" "<<cornerNr+2*nrCellsY+1<<" "<<cornerNr+3*nrCellsY+1<<" "<<cornerNr+2*nrCellsY<<" "<<cornerNr+nrCellsY<<endl
					<<"  nrsprings "<<6
					<<" numbers "
					<<wall[2*j+1][2*i]<<" "<<wall[2*j+2][2*i]<<" "<<wall[2*j+2][2*i+1]<<" "<<wall[2*j+2][2*i+2]<<" "<<wall[2*j+1][2*i+2]<<" "<<wall[2*j+1][2*i+1]<<endl;
				cornerNr++;
			}
			cornerNr += nrCellsY+1;
			wallNr += nrCellsY;
		}
	}
	outFile<<"nrcells "<<cellNr<<endl;
	outFile<<"-end-configuration-"<<endl<<endl;

	
	outFile<<"-northface-"<<endl;
	outFile<<"nrcorners "<<nrCellsX+1<<endl;
	outFile<<"numbers ";
	for(i=2*nrCellsY+1; i<nrCorners; )
	{
		outFile<<i<<" ";	i+=nrCellsY+1;
		if(i<=nrCorners) {outFile<<i<<" ";	i+=3*nrCellsY+1;}
	}
	outFile<<endl<<"-end-northface-"<<endl<<endl;

	
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrCellsY<<endl;
	outFile<<"numbers ";
	for(i=nrCorners; i>nrCorners-nrCellsY; --i)
	{
		outFile<<i<<" ";
	}
	outFile<<endl<<"-end-eastface-"<<endl<<endl;


	outFile<<"-southface-"<<endl;
	outFile<<"nrcorners "<<nrCellsX+1<<endl;
	outFile<<"numbers ";
	for(i=nrCellsY+1; i<nrCorners; )
	{
		outFile<<i<<" ";	i+=nrCellsY+1;
		if(i<=nrCorners) {outFile<<i<<" ";	i+=3*nrCellsY+1;}
	}
	outFile<<endl<<"-end-southface-"<<endl<<endl;


	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrCellsY<<endl;
	outFile<<"numbers ";
	for(i=1; i<nrCellsY+1; ++i)
	{
		outFile<<i<<" ";
	}
	outFile<<endl<<"-end-westface-"<<endl<<endl;

		
	outFile<<"-crosssection-"<<endl;
	int middle = (int)floor(nrCellsX/2.0)+1;
	int firstWall, firstSecWall;
	int nrCrossSprings;
	int nrCrossSecSprings;
	if(middle%2==0) // small column
	{
		nrCrossSprings = nrCellsY;
		nrCrossSecSprings = 2*(nrCellsY - 1);
		int half = floor(middle/2)-1;
		firstWall = half*(4*nrCellsY + nrCellsY+1 + nrCellsY) + 4*nrCellsY + nrCellsY + 2;
		firstSecWall = half*(2*nrCellsY-1)*2 + nrCellsY*2 + 1;
	}else	// big column
	{
		nrCrossSprings = nrCellsY+1;
		nrCrossSecSprings = 4*nrCellsY;
		int half = floor(middle/2);
		firstWall = half*(4*nrCellsY + nrCellsY+1 + nrCellsY) + 1;
		firstSecWall = half*(2*nrCellsY-1)*2 + 1;
	}
	outFile<<"nrsprings "<<nrCrossSprings<<endl;
	outFile<<"numbers ";
	for(i=0; i<nrCrossSprings; ++i)
	{
		outFile<<firstWall++<<" ";
	}
	outFile<<endl;
	outFile<<"nrsecundarysprings "<<nrCrossSecSprings<<endl;
	outFile<<"secundarynumbers ";
	for(i=0; i<nrCrossSecSprings;)
	{
/*		if(i++<nrCrossSecSprings){outFile<<firstSecWall++<<" ";}
		if(i++<nrCrossSecSprings){outFile<<firstSecWall++<<" ";}
		firstSecWall += 2;
		if(i++<nrCrossSecSprings){outFile<<firstSecWall++<<" ";}
		if(i++<nrCrossSecSprings){outFile<<firstSecWall++<<" ";}
		*/
		if(i++<nrCrossSecSprings){outFile<<firstSecWall++<<" ";}
	}
	outFile<<endl<<"-end-crosssection-"<<endl<<endl;

	outFile.close();
}

///////////////////////////////////////////////////////////
// o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o //
///////////////////////////////////////////////////////////

void transversal()
{
		const char* nameOut = "data/inputFileT.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);

	// switch x and y for the imput
	double cellDimX = cllDimX;
	double cellDimY = cllDimY;
	const int		nrCellsX = nrCllsXT;
	const int		nrCellsY = nrCllsYT;

	int dxWidth = 6;
	double dx = cellDimX/(double)(dxWidth+2); 
	double dy = cellDimY/(double)(2);
	int		nrDX = (dxWidth+1)*nrCellsX+1;
	int		nrDY = 2*nrCellsY;

	int		bigColumnDY = nrDY+1;
	int		smallColumnDY = nrDY;
	int		columnCode = 1;

	int cornerNr = 0;
	double x,y;
	outFile<<"-cornerpoints-"<<endl;
	for(int ix=0; ix<=nrDX; )
	{
		x = ix*dx;
		if((columnCode==1)||(columnCode==4))
		{
			for(int iy=1; iy<=smallColumnDY; iy += 2)
			{
				 y = iy*dy;
				 outFile<<++cornerNr<<" "<<y<<" "<<x<<endl;
			}
			if(columnCode==1){columnCode=2; ix+=1;}
			else{/*(columnCode==4)*/columnCode=1; ix+=dxWidth;}
		} else { // (columnCode==2)||(columnCode==3)
			for(int iy=0; iy<=bigColumnDY; iy += 2)
			{
				 y = iy*dy;
				 outFile<<++cornerNr<<" "<<y<<" "<<x<<endl;
			}
			if(columnCode==2){columnCode=3; ix+=dxWidth;}
			else{/*(columnCode==3)*/columnCode=4; ix+=1;}
		}
	}
	int nrCorners = cornerNr;
	outFile<<"nrcorners "<<nrCorners<<endl;
	outFile<<"-end-cornerpoints-"<<endl<<endl;

////////////////////////////////////////////////////////////////////
// get the springs right
////////////////////////////////////////////////////////////////////
	int springNr = 0;
	columnCode = 1;
	int rightBegin = 1;
	int left;
	int	right;
	outFile<<"-springs-"<<endl;
	int i;
	for(ix=0; ix<nrDX;  )
	{
		if(columnCode==1)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY;
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<right++<<" "<<left  <<endl;
				 outFile<<++springNr<<" "<<right  <<" "<<left++<<endl;
			}
			{columnCode=2; ix+=1;}
		}else{
		if(columnCode==2)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY+1;
			right = rightBegin;
			for(i=0; i<nrCellsY+1; ++i)
			{
				 outFile<<++springNr<<" "<<right++<<" "<<left++<<endl;
			}
			{columnCode=3; ix+=dxWidth;}
		}else{
		if(columnCode==3)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY+1;
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<right  <<" "<<left++<<endl;
				 outFile<<++springNr<<" "<<right++<<" "<<left  <<endl;
			}
			{columnCode=4; ix+=1;}
		}else{
		if(columnCode==4)
		{
			left = rightBegin;
			rightBegin = left + nrCellsY; // one less than the others!
			right = rightBegin;
			for(i=0; i<nrCellsY; ++i)
			{
				 outFile<<++springNr<<" "<<right++<<" "<<left++<<endl;
			}
			{columnCode=1; ix+=dxWidth;}
		}
		} // end if(columnCode==1)else
		} // end if(columnCode==2)else
		} // end if(columnCode==3)else
	}
	int nrSprings = springNr;
	outFile<<"nrsprings "<<nrSprings<<endl;
	outFile<<"-end-springs-"<<endl<<endl;

	////////////////////////////////////////////////////////////////////
// get the secundary springs right
////////////////////////////////////////////////////////////////////
	springNr = 0;
	int wi, w, sw, nw, se, ne, e;
	bool thickColumn = true;
	wi = w = 1;
	sw = w + nrCellsY;
	nw = sw + 1;
	se = sw + nrCellsY + 1;
	ne = se + 1;
	e = se + nrCellsY + 1;

	outFile<<"-secundarysprings-"<<endl;
	int cellColumnX;
	int cellRowY = 0;
	for(cellColumnX=0; cellColumnX<nrCellsX;  )
	{
		while(cellRowY<nrCellsY)
		{
//			outFile<<++springNr<<" "<<sw<<" "<<e <<endl;
//			outFile<<++springNr<<" "<<se<<" "<<w <<endl;
			outFile<<++springNr<<" "<<sw<<" "<<nw<<endl;
			outFile<<++springNr<<" "<<se<<" "<<ne<<endl;
//			outFile<<++springNr<<" "<<w <<" "<<ne<<endl;
//			outFile<<++springNr<<" "<<e <<" "<<nw<<endl;

			sw = nw;
			se = ne;
			w++;
			e++;
			nw++;
			ne++;
			cellRowY++;
			
		}
		

		if(thickColumn){
			wi += 2*nrCellsY + 2;
			w = wi;
			sw = w + nrCellsY;
			nw = sw + 1;
			se = sw + nrCellsY;
			ne = se + 1;
			e  = se + nrCellsY + 1;
			cellRowY = 1;
		} else {
			wi += 2*nrCellsY;
			w = wi;
			sw = w + nrCellsY;
			nw = sw + 1;
			se = sw + nrCellsY + 1;
			ne = se + 1;
			e  = se + nrCellsY + 1;
			cellRowY = 0;
		}
		cellColumnX++;
		thickColumn = !thickColumn;
	}
	int nrSecundarySprings = springNr;
	outFile<<"nrsprings "<<nrSecundarySprings<<endl;
	outFile<<"-end-secundarysprings-"<<endl<<endl;
	

////////////////////////////////////////////////////////////////////
// get the cells right
////////////////////////////////////////////////////////////////////

	int wallNr = 1;
	int wallColumn=1;
	int wall[2*nrCellsY][2*nrCellsX+1];
	for(int j=0; j<2*nrCellsX+1; ++j){
		if((wallColumn==1)||(wallColumn==3))
		{
			for(i=0; i<2*nrCellsY; ++i)
			{
				wall[i][j] = wallNr++;
			}
			wallColumn++;
		}else
		{
			if(wallColumn==2)
			{
				for(i=0; i<2*nrCellsY; ++i)
				{
					if(i%2!=0)
					{
						wall[i][j] = wallNr;
					}else
					{
						wall[i][j] = wallNr++;
					}	
				}
				wallNr++;
				wallColumn = 3;
			}else
			{
				if(wallColumn==4)
				{
					wall[0][j] = 0;
					for(i=1; i<2*nrCellsY-1; ++i)
					{
						if(i%2==0)
						{
							wall[i][j] = wallNr;
						}else
						{
							wall[i][j] = wallNr++;
						}	
					}
					wallNr++;
					wall[2*nrCellsY-1][j] = 0;
					wallColumn = 1;
				}
			}
		}
	}

	outFile<<"-configuration-"<<endl;
	// first row of cells
	int cellNr=0;
	cornerNr = 1;

	for(i=0; i<nrCellsX; ++i)
	{
		if(i%2==0) // a big column of nrCellsY
		{
			for(int j=0; j<nrCellsY; ++j)
			{
				outFile<<++cellNr
					<<" nrcorners "<<6
					<<" numbers "
					<<cornerNr<<" "<<cornerNr+nrCellsY+1<<" "<<cornerNr+2*nrCellsY+2<<" "<<cornerNr+3*nrCellsY+2<<" "<<cornerNr+2*nrCellsY+1<<" "<<cornerNr+nrCellsY<<endl
					<<"  nrsprings "<<6
					<<" numbers "
					<<wall[2*j][2*i]<<" "<<wall[2*j+1][2*i]<<" "<<wall[2*j+1][2*i+1]<<" "<<wall[2*j+1][2*i+2]<<" "<<wall[2*j][2*i+2]<<" "<<wall[2*j][2*i+1]<<endl;
				cornerNr++;
			}
			cornerNr += nrCellsY+2;
			wallNr += nrCellsY+2;
		}else{ // a small column of nrCellsY-1
			for(int j=0; j<nrCellsY-1; ++j)
			{
				outFile<<++cellNr
					<<" nrcorners "<<6
					<<" numbers "
					<<cornerNr<<" "<<cornerNr+nrCellsY+1<<" "<<cornerNr+2*nrCellsY+1<<" "<<cornerNr+3*nrCellsY+1<<" "<<cornerNr+2*nrCellsY<<" "<<cornerNr+nrCellsY<<endl
					<<"  nrsprings "<<6
					<<" numbers "
					<<wall[2*j+1][2*i]<<" "<<wall[2*j+2][2*i]<<" "<<wall[2*j+2][2*i+1]<<" "<<wall[2*j+2][2*i+2]<<" "<<wall[2*j+1][2*i+2]<<" "<<wall[2*j+1][2*i+1]<<endl;
				cornerNr++;
			}
			cornerNr += nrCellsY+1;
			wallNr += nrCellsY;
		}
	}
	outFile<<"nrcells "<<cellNr<<endl;
	outFile<<"-end-configuration-"<<endl<<endl;

	
	outFile<<"-eastface-"<<endl;
	outFile<<"nrcorners "<<nrCellsX+1<<endl;
	outFile<<"numbers ";
	for(i=2*nrCellsY+1; i<nrCorners; )
	{
		outFile<<i<<" ";	i+=nrCellsY+1;
		if(i<=nrCorners) {outFile<<i<<" ";	i+=3*nrCellsY+1;}
	}
	outFile<<endl<<"-end-eastface-"<<endl<<endl;

	
	outFile<<"-northface-"<<endl;
	outFile<<"nrcorners "<<nrCellsY<<endl;
	outFile<<"numbers ";
	for(i=nrCorners; i>nrCorners-nrCellsY; --i)
	{
		outFile<<i<<" ";
	}
	outFile<<endl<<"-end-northface-"<<endl<<endl;


	outFile<<"-westface-"<<endl;
	outFile<<"nrcorners "<<nrCellsX+1<<endl;
	outFile<<"numbers ";
	for(i=nrCellsY+1; i<nrCorners; )
	{
		outFile<<i<<" ";	i+=nrCellsY+1;
		if(i<=nrCorners) {outFile<<i<<" ";	i+=3*nrCellsY+1;}
	}
	outFile<<endl<<"-end-westface-"<<endl<<endl;


	outFile<<"-southface-"<<endl;
	outFile<<"nrcorners "<<nrCellsY<<endl;
	outFile<<"numbers ";
	for(i=1; i<nrCellsY+1; ++i)
	{
		outFile<<i<<" ";
	}
	outFile<<endl<<"-end-southface-"<<endl<<endl;

		
	outFile<<"-crosssection-"<<endl;
	int firstWall = nrCellsY;
	
	outFile<<"nrsprings "<<nrCellsX + 1<<endl;
	outFile<<"numbers ";
	i=firstWall;
	bool bigColumn = true;
	outFile<<i<<" ";
	for(; i<nrSprings;)
	{
		if(bigColumn){i += 3*nrCellsY+1;}
		else{i += 3*nrCellsY;}
		if(i<nrSprings){outFile<<i<<" ";}
		bigColumn = !bigColumn;
	}
	outFile<<endl;
	int halfNrCellsY = (int)floor(nrCellsY/2.0);
	bool belowHalf;
	bigColumn = false;
	if(nrCellsY%2==0){ firstWall = halfNrCellsY * 2 -1; belowHalf = true;}
	else { firstWall = halfNrCellsY*2 + 1; belowHalf = false; }
	outFile<<"nrsecundarysprings "<<2*nrCellsX <<endl;
	/*	if(nrCellsY%2==0){ firstWall = halfNrCellsY * 6 -1; belowHalf = true;}
	else { firstWall = halfNrCellsY*6 + 1; belowHalf = false; }
	outFile<<"nrsecundarysprings "<<4*nrCellsX <<endl;
	*/
	outFile<<"secundarynumbers ";
	for(i=firstWall; i+2<nrSecundarySprings;)
	{
		 outFile<<i++<<" ";//outFile<<i+1<<" "; outFile<<i+2<<" ";
		if(i<nrSecundarySprings){outFile<<i<<" ";}
		bigColumn = !bigColumn;
		if(bigColumn){
			if(belowHalf){ i += halfNrCellsY*2 + (halfNrCellsY-1)*2 + 1; belowHalf=false;}
			else{ i += (halfNrCellsY*2 ) + (halfNrCellsY)*2 - 1; belowHalf=true; }
		}
		else{
			if(belowHalf){ i += halfNrCellsY*2 + (halfNrCellsY)*2 + 1; belowHalf=false;}
			else{ i += ((halfNrCellsY-1)*2) + (halfNrCellsY)*2 - 1; belowHalf=true; }
		}
		
		if(i<nrSecundarySprings){/*outFile<<i-1<<" "; outFile<<i-2<<" ";*/outFile<<i++<<" ";}
		if(i<nrSecundarySprings){outFile<<i<<" ";}
		bigColumn = !bigColumn;
		if(bigColumn){
			if(belowHalf){ i += halfNrCellsY*2 + (halfNrCellsY-1)*2 + 1; belowHalf=false;}
			else{ i += (halfNrCellsY*2 ) + (halfNrCellsY)*2 - 1; belowHalf=true; }
		}
		else{
			if(belowHalf){ i += halfNrCellsY*2 + (halfNrCellsY)*2 + 1; belowHalf=false;}
			else{ i += ((halfNrCellsY-1)*2) + (halfNrCellsY)*2 - 1; belowHalf=true; }
		}
	}
		/*if(bigColumn){
			if(belowHalf){ i += halfNrCellsY*6 + (halfNrCellsY+1)*6 + 1; belowHalf=false;}
			else{ i += (halfNrCellsY*6 + 4) + (halfNrCellsY)*6 - 1; belowHalf=true; }
		}
		else{
			if(belowHalf){ i += halfNrCellsY*6 + (halfNrCellsY)*6 + 1; belowHalf=false;}
			else{ i += ((halfNrCellsY-1)*6 + 4) + (halfNrCellsY)*6 - 1; belowHalf=true; }
		}*/

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