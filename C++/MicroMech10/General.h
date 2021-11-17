// The debugger can't handle symbols more than 255 characters long.
// STL often creates symbols longer than that.
// When symbols are longer than 255 characters, the warning is disabled.
// #pragma warning(disable:4786)

//////////////////////////////////////////////////////////////////
#ifndef GENERAL_H												//
#define GENERAL_H												//
// begin of General.h											//
//////////////////////////////////////////////////////////////////

#include <math.h>
#include <iostream.h> 

// Pi
#define PI 3.141592653589793238462643

// gravity acceleration in the negative z-direction
#define GRAVITY 9.81
					

class CTriple
{
public:
	double	x,	// x-value of the CTriple
			y,  // y-value of the CTriple
			z;	// z-value of the CTriple
	CTriple(){}
	CTriple(double newx, double newy, double newz)
			:x(newx),   y(newy),   z(newz){}

	CTriple operator+ (CTriple &B)
	{ // sum of two vectors
		return CTriple(x + B.x, y + B.y, z + B.z); }
	CTriple operator- (CTriple &B)
	{ // difference of two vectors	
		return CTriple(x - B.x, y - B.y, z - B.z); }
	double operator* (CTriple &B)
	{ // dot product between two vectors
		return x * B.x + y * B.y + z * B.z;}
	CTriple operator* (double b)
	{ // vector multiplied by a scalar
		return CTriple(b * x, b * y, b * z);}
	friend CTriple operator* (double b, CTriple &A)
	{ // scalar multiplied by a vector
		return CTriple(b * A.x, b * A.y, b * A.z);	}
	CTriple operator/ (double b)
	{ // vector divided by a scalar
		if (b != 0)return CTriple(x / b, y / b, z / b);	
		else {cout<<"error: division by zero"; return CTriple(0,0,0);}}
};															

#define ZEROTRIPLE CTriple(0,0,0)

static int Dir(CTriple &A)
	{// measure to detect the switching of direction
		if (A.x+A.y+A.z<0) return -1;
		else	 return 1;}

static double Distance(CTriple &A, CTriple &B)
{ // distance between two CTriples
	return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y)+(A.z-B.z)*(A.z-B.z));
}

static CTriple VectorProd(CTriple &A, CTriple &B)
{ 
	return CTriple(A.y*B.z-A.z*B.y,A.z*B.x-A.x*B.z,A.x*B.y-A.y*B.x);
}

static double VectorLength(CTriple &A)
{
	return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);
}
	
static CTriple NormVector(CTriple &A, CTriple &B)
{
	CTriple Q(VectorProd(A,B));
	double IIQII	= VectorLength(Q);
	if (IIQII==0)		return CTriple(0,0,0);
	else Q = Q/IIQII;	return Q;
}

class CPair
{
public:
	double	x,	// x-value of the CPair
			y;  // y-value of the CPair
	CPair(){}
	CPair(double newx, double newy)
			:x(newx),   y(newy){}
	CPair(CPair & newP)
		:x(newP.x),	y(newP.y){}

	CPair  operator- ()	const
	{ // unary minus
		return CPair(-x, -y);			}
	CPair operator+ (CPair &B)
	{ // sum of two vectors
		return CPair(x + B.x, y + B.y); }
	CPair operator- (CPair &B)
	{ // difference of two vectors	
		return CPair(x - B.x, y - B.y); }
	double operator* (CPair &B)
	{ // dot product between two vectors
		return x * B.x + y * B.y;		}
	CPair operator* (double b)
	{ // vector multiplied by a scalar
		return CPair(b * x, b * y);		}
	friend CPair operator* (double b, CPair &A)
	{ // scalar multiplied by a vector
		return CPair(b * A.x, b * A.y);	}
	CPair operator/ (double b)
	{ // vector divided by a scalar
		if (b != 0)return CPair(x / b, y / b);	
		else {cout<<"error: division by zero"; return CPair(0,0);}}
};															

#define ZEROPAIR CPair(0,0)

static int Dir(CPair &A)
	{// measure to detect the switching of direction
		if (A.x+A.y < 0) return -1;
		else	 return 1;}

static double Distance(CPair &A, CPair &B)
{ // distance between two CPairs
	return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y));
}

static CTriple VectorProd(CPair &A, CPair &B)
{ 
	return CTriple(0,0,A.x*B.y-A.y*B.x);
}

static double VectorLength(CPair &A)
{
	return sqrt(A.x*A.x+A.y*A.y);
}
	
static CTriple NormVector(CPair &A, CPair &B)
{
	CTriple Q(VectorProd(A,B));
	double IIQII	= VectorLength(Q);
	if (IIQII==0)		return CTriple(0,0,0);
	else Q = Q/IIQII;	return Q;
}

static double Angle(CPair &A)
{
	if(A.x>0)
	{
		if(A.y>0)
			return atan(A.y/A.x);
		else
			return PI-atan(A.y/A.x);
	}
	else // A.x<0
	{
		if(A.y>0)
			return atan(A.y/A.x);
		else
			return atan(A.y/A.x)-PI;
	}
}

static double Angle(CPair &A, CPair &B)
{
	// the corner is described by the points acb, where c is the center point
	// a is CW from c and b is CCW from c, alfa is the angle CCW is positive
	/*	
     *         _ b
     *         /| 
     *    B  /
     *     / alfa
     *   /_)________\ a
     * c      A     /  
	 *
	 *                (A dot B)     
	 *   cos(alfa) =  ---------
	 *                (|A||B|)
	 */
	// acos: arc cosinus returns the angle in radians (while openGL works with degrees)
	return Angle(B)-Angle(A);
}

class CQuater
{
public:
	double	s,	// s-value of the quaternion
			i,  // i-value of the quaternion
			j,	// j-value of the quaternion
			k;	// k-value of the quaternion
	CQuater(){}
	CQuater(double news, double newi, double newj, double newk)
			:s(news), i(newi), j(newj), k(newk){}

	void Normalize() {double length = sqrt(s*s+i*i+j*j+k*k); 
					  s=s/length; i=i/length; j=j/length; k=k/length; }
	CQuater operator+ (CQuater &B)
	{ // sum of two quaternions
		return CQuater(s + B.s, i + B.i, j + B.j, k + B.k); }
	CQuater operator* (CQuater &B)
	{ // product between two quaternions
		return CQuater(s*B.s-i*B.i-j*B.j-k*B.k,
					   i*B.s+s*B.i-k*B.j+j*B.k,
					   j*B.s+k*B.i-s*B.j+i*B.k,
					   k*B.s+j*B.i-i*B.j+s*B.k);}
	friend CQuater operator* (CTriple &A, CQuater &B)
	{ // product between a triple and a quaternion
		return CQuater(-A.x*B.i-A.y*B.j-A.z*B.k,
					   A.x*B.s-A.z*B.j+A.y*B.k,
					   A.y*B.s+A.z*B.i+A.x*B.k,
					   A.z*B.s+A.y*B.i-A.x*B.j);}
	CQuater operator* (double b)
	{ // quaternion multiplied by a scalar
		return CQuater(b * s, b * i, b * j, b * k);}
	friend CQuater operator* (double b, CQuater &A)
	{ // scalar multiplied by a quaternion
		return CQuater(b * A.s, b * A.i, b * A.j, b * A.k);	}
};		

class CMatrix
{
public:
	double	xx,	xy,	xz,
			yx,	yy,	yz,
			zx,	zy,	zz;
	CMatrix(){}
	CMatrix(double newxx, double newxy, double newxz, 
			double newyx, double newyy, double newyz,
			double newzx, double newzy, double newzz)
			:xx(newxx),   xy(newxy),	xz(newxz),
			 yx(newyx),   yy(newyy),	yz(newyz),
			 zx(newzx),   zy(newzy),	zz(newzz){}

	CMatrix operator+ (CMatrix &B)
	{ // sum of two matrices
		return CMatrix(xx + B.xx, xy + B.xy, xz + B.xz,
					   yx + B.yx, yy + B.yy, yz + B.yz,
					   zx + B.zx, zy + B.zy, zz + B.zz); }
	CMatrix operator- (CMatrix &B)
	{ // difference of two matrices	
		return CMatrix(xx - B.xx, xy - B.xy, xz - B.xz,
					   yx - B.yx, yy - B.yy, yz - B.yz,
					   zx - B.zx, zy - B.zy, zz - B.zz); }
	CMatrix operator* (CMatrix &B)
	{ // product between two matrices
		return CMatrix(xx*B.xx+xy*B.yx+xz*B.zx, xx*B.xy+xy*B.yy+xz*B.zy, xx*B.xz+xy*B.yz+xz*B.zz,
					   yx*B.xx+yy*B.yx+yz*B.zx, yx*B.xy+yy*B.yy+yz*B.zy, yx*B.xz+yy*B.yz+yz*B.zz,
					   zx*B.xx+zy*B.yx+zz*B.zx, zx*B.xy+zy*B.yy+zz*B.zy, zx*B.xz+zy*B.yz+zz*B.zz); }
	CTriple operator* (CTriple &B)
	{ // matrix multiplied by a vector
		return CTriple(xx*B.x+xy*B.y+xz*B.z, yx*B.x+yy*B.y+yz*B.z, zx*B.x+zy*B.y+zz*B.z); }
	friend CTriple operator* (CTriple &B, CMatrix &A)
	{ // vector multiplied by a matrix
		return CTriple(A.xx*B.x+A.yx*B.y+A.zx*B.z, 
					   A.xy*B.x+A.yy*B.y+A.zy*B.z, 
					   A.xz*B.x+A.yz*B.y+A.zz*B.z); }
	CMatrix operator* (double b)
	{ // matrix multiplied by a scalar
		return CMatrix(b*xx,	b*xy,	b*xz,
					   b*yx,	b*yy,	b*yz,
					   b*zx,	b*zy,	b*zz); }
	friend CMatrix operator* (double b, CMatrix A)
	{ // scalar multiplied by a matrix
		return CMatrix(b*A.xx,	b*A.xy,	b*A.xz,
					   b*A.yx,	b*A.yy,	b*A.yz,
					   b*A.zx,	b*A.zy,	b*A.zz); }
	CMatrix operator/ (double b)
	{ // matrix divided by a scalar
		if (b != 0)return CMatrix(xx/b, xy/b, xz/b,
								  yx/b, yy/b, yz/b,
								  zx/b, zy/b, zz/b);	
		else {cout<<"error: division by zero"; return CMatrix(0,0,0,0,0,0,0,0,0);}}
};	

#define ZEROMATRIX CMatrix(0,0,0,0,0,0,0,0,0)
#define UNITMATRIX CMatrix(1,0,0,0,1,0,0,0,1)

static CMatrix Transpose(CMatrix &M)
{ // transpose of a matrix
	return CMatrix(M.xx,	M.yx,	M.zx,
				   M.xy,	M.yy,	M.zy,
				   M.xz,	M.yz,	M.zz); 
}

static CMatrix Inverse(CMatrix &A)
{	// C is the coefficient matrix of A
	if ((A.xy==0)&&(A.xz==0)&&(A.yx==0)&&(A.yz==0)&&(A.zx==0)&&(A.zy==0)
		&&(A.xx!=0)&&(A.yy!=0)&&(A.zz!=0))
	{
		return CMatrix(1/A.xx, 0, 0, 0, 1/A.yy, 0, 0, 0, 1/A.zz);
	}
	else
	{
		CMatrix C(A.yy*A.zz-A.yz*A.zy, -A.yx*A.zz+A.yz*A.zx, A.yx*A.zy-A.yy*A.zx,
				  -A.xy*A.zz+A.xz*A.zy, A.xx*A.zz-A.xz*A.zx, -A.xx*A.zy+A.xy*A.zx,
				  A.xy*A.yz-A.xz*A.yy, -A.xx*A.yz+A.xz*A.yx, A.xx*A.yy-A.xy*A.yx);
		double IAI	= A.xx*C.xx+A.xy*C.xy+A.xz*C.xz;
		if (IAI==0)		return ZEROMATRIX;
		else return CMatrix(C.xx/IAI, C.yx/IAI, C.zx/IAI,
							C.xy/IAI, C.yy/IAI, C.zy/IAI,
							C.xz/IAI, C.yz/IAI, C.zz/IAI);
	}
}	

static CMatrix Quater_to_Matrix(CQuater &Q)
{
	Q.Normalize();
	return CMatrix(1-2*Q.j*Q.j-2*Q.k*Q.k,	2*Q.i*Q.j-2*Q.s*Q.k,	2*Q.i*Q.k+2*Q.s*Q.j,
				   2*Q.i*Q.j+2*Q.s*Q.k,		1-2*Q.i*Q.i-2*Q.k*Q.k,	2*Q.j*Q.k-2*Q.s*Q.i,
				   2*Q.i*Q.k-2*Q.s*Q.j,		2*Q.j*Q.k+2*Q.s*Q.i,	1-2*Q.i*Q.i-2*Q.j*Q.j);
}

static CQuater Matrix_to_Quater(const CMatrix &M)
{
	CQuater	Q;
	double	tr, s;
	
	tr = M.xx+M.yy+M.zz;

	if (tr>=0)
	{
		s = sqrt(tr+1);
		Q.s = 0.5*s;
		s = 0.5/s;
		Q.i = (M.zy-M.yz)*s;
		Q.j = (M.xz-M.zx)*s;
		Q.k = (M.yx-M.xy)*s;
	}
	else
	{
		int i=0;
		if (M.yy > M.xx) {i=1; if (M.zz > M.yy) i=2;}
		else if (M.zz > M.xx) i=2;
		switch (i)
		{
		case 0:
			s = sqrt((M.xx-(M.yy+M.zz))+1);
			Q.i = 0.5*s;
			s = 0.5/s;
			Q.j = (M.xy-M.yx)*s;
			Q.k = (M.zx-M.xz)*s;
			Q.s = (M.zy-M.yz)*s;
			break;

		case 1:
			s = sqrt((M.yy-(M.zz+M.xx))+1);
			Q.j = 0.5*s;
			s = 0.5/s;
			Q.k = (M.yz-M.zy)*s;
			Q.i = (M.xy-M.yx)*s;
			Q.s = (M.xz-M.zx)*s;
			break;

		case 2:
			s = sqrt((M.zz-(M.xx+M.yy))+1);
			Q.k = 0.5*s;
			s = 0.5/s;
			Q.i = (M.zx-M.xz)*s;
			Q.j = (M.yz-M.zy)*s;
			Q.s = (M.yx-M.xy)*s;
			break;
		}
	}
	return Q;
}

//////////////////////////////////////////////////////////////////
#endif															//
// end of General.h												//
//////////////////////////////////////////////////////////////////
