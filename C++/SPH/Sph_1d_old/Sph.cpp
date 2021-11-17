/****************************************************************************
 *                                                                            
 *                                                                            
 *                                                                            
 *                                 SPH_1D                                      
 *                                                                            
 *                               KEVIN OLSON                                  
 *                           USRA and NCCS, GSFC                              
 *                                                                            
 *                              April, 1993                                    
 *                                                                            
 *     Smooth Particle Hydrodynamics is a gridless particle method for        
 *     solving the Lagrangian equations of hydrodynamics.  This               
 *     implementation includes variable smoothing lengths using both the      
 *     algorithm described in Benz (1990) and that described in Steinmetz     
 *     and Muller (1992).  Also included is a routine to compute an           
 *     additional thermal diffusion term which is added to the energy         
 *     equation (Monaghan 1988).  I have found that in problems where         
 *     the initial conditions are discontinuous, large deviations in the      
 *     energy profile occur at the contact discontinuity which forms near     
 *     the initial discontinuity.  Addition of a thermal diffusion term       
 *     removes this problem for the double shock tube problem (Woodward and   
 *     Colella, 1984).  Also included in the model are routines for           
 *     both the molecular and bulk viscosities.                               
 *                                                                            
 *     This implementation of SPH has been written in Maspar Fortran to run   
 *     on the Maspar MP-1 at GSFC.                                            
 *                                                                             
 *     The following references where used extensively in building this       
 *     code:                                                                  
 *                                                                            
 *      Benz, W., 1990, in Numerical Modelling of Nonlinear Stellar           
 *         Pulsations, Problems and Prospects, ed. J.R. Buchler, Kulwer       
 *         Acdemics Publishers, Dordecht, Boston, London                      
 *                                                                            
 *      Hernquist, L., and Katz, N., 1989, ApJ. Suppl., Vol. 70, 419          
 *                                                                           
 *      Monaghan, J.J., 1985, Comp. Phys. Rep., Vol. 3, 71                    
 *                                                                            
 *      Monaghan, J.J., 1988, "SPH Meets the Shocks of Noh", unpublished      
 *         preprint                                                           
 *                                                                            
 *      Steinmetz, M., and Muller, E., 1992, preprint sub. to A&A             
 *                                                                            
 *      Woodward, P., and Colella, P., 1984, J. Comp. Phys., Vol. 54, 115     
 *                                                                            
 *                                                                            
 ***************************************************************************** 
 */

#include <math.h>
#include <iostream>
#include <fstream>

#include "Noh_int.h"
#include "R4_sort.h"
#include "Sod_const.h"
#include "Sod_int.h"

#if defined(_WIN32)
#include <sys/timeb.h>
#define CLK_TCK 1000
#else
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#endif

float TimeDelta(void)
{
/* 
	timedelta: returns the number of seconds that have elapsed since
	the previous call to the function. 
*/
    static long begin = 0;
    static long finish, difference;
#if defined(_WIN32)
    static struct timeb tb;
    ftime(&tb);
    finish = tb.time*1000+tb.millitm;
#else
    static struct tms tb;
    finish = times(&tb);
#endif
    difference = finish - begin;
    begin = finish;
    return (float)difference/(float)CLK_TCK;
}

using namespace std;

	// set nsize to the total no. of particles plus 100 for the boundries
	const int nsize = 400 + 100; 

	double	xo[nsize], vxo[nsize], uo[nsize];
	double	vxt[nsize];
	double	deno[nsize];
	double	alpha, beta;
	double	c1, c2, tend;
	double	alphab, betab, g1, g2;
	double	hfac;
	double	ho[nsize];
	double	sound[nsize];
	double	maxmu[nsize];
	double	tote, toteo;

	int		isort[nsize];
            
	/******************************************************************************
	 *                                                                           
	 *     VARIABLE LIST                                               
	 */
	
	//     xx,vx,ax:  -> position,velocity, and acceleration of a particle       
	double xx[nsize], vx[nsize], ax[nsize];

	//     u,au:      -> energy per unit mass and its rate of change, du/dt      
	double u[nsize], au[nsize];
                                                                          
	//     den,p,mass:-> density, pressure, and mass of a particle            
	double den[nsize], p[nsize], mass[nsize];

	//     hh,mass:   -> the smoothing length and mass of a particle             
	double hh[nsize];

	//     t,delt:    -> the time and the time step
	double t,delt;

	//     divv:      -> the divergence of v
	double divv[nsize];

	//     thermdiff: -> the additional thermal diffusion term added to du/dt    
	double thermdiff[nsize];

	//     gamma:     -> the ratio of specific heats                             
	double gamma;

	//     cour:      -> the courant number                                      
	double cour;

	//     partno:    -> the number of a particle in the list                    
	int partno[nsize];

	//	   xmin, xmax:-> the minimum and maximum values of x, 
	//					 they define the domain of the problem 
	double xmin, xmax;

	//	   n:		  -> the total number of particles 
	const int n = nsize;

	//	   nsteps:	  -> a counter which counts the number of time steps taken
	int nsteps;

	// used for the sorting operations
	double vxTemp[nsize];
	double axTemp[nsize];
	double pTemp[nsize];
	double denTemp[nsize];
	double hhTemp[nsize];
	double massTemp[nsize];
	double uTemp[nsize];
	double auTemp[nsize];
	double xoTemp[nsize];
	double vxoTemp[nsize];
	double uoTemp[nsize];
	double hoTemp[nsize];

	// other arrays used to hold temporary data in the subroutines
	double w[nsize];
	double w2[nsize];
	double gradw2[nsize];
	double q[nsize];
	double xt[nsize];
	double axt[nsize];
	double aut[nsize];
	double ht[nsize];
	double delti[nsize];
	double pt[nsize];
	double ut[nsize];
	double masst[nsize];
	double dent[nsize];
	double denst[nsize];
	double divvt[nsize];
	double soundt[nsize];
	double maxmut[nsize];
	double r[nsize];
	double v[nsize];
	double h[nsize];
	double gradw[nsize];
	double qi[nsize];
	double qj[nsize];
	double cij[nsize];
	double denij[nsize];
	double artv[nsize];
	double eta[nsize];
	double mu[nsize];
	double pforce[nsize];
	double thermdifft[nsize];
	double uu[nsize];
	double qqi[nsize];
	double qqj[nsize];
	double old_u[nsize];
	double old_ut[nsize];
	double dsmt[nsize];
	double densmooth[nsize];
	double dsmoothdt[nsize];
	double avh[nsize];

	/*     The following are switches to be set to determine if and which artificial
	 *     viscosity is to be used, as well as if the thermal diffusion term in the 
	 *     energy equation is to be included.  A value of 1 indicates that the 
	 *     switch is on.
	 */
	bool  ibulk=0;          // the bulk viscosity switch

	bool  inorm=1;          // the molecular viscosity switch
      
	bool  itherm=0;         // the thermal diffusion switch
							// note: do not use thermal diffusion with entropy

	bool  ientro=0;         // specifies whether or not the entropic equation
							// is to be used in place of the energy equation

	bool  ientro2=0;        // specifies whether entropy is used in lieu of the
							// of the force equation

	bool  irefl=1;          // specifies whether reflecting boundries are used

	/*     the following, twr1-twr8, define the times at which output will 
	 *     be produced
	 */
	double twr1 = 0.00001;
	double twr2 = 0.000016;
	double twr3 = 0.000026;
	double twr4 = 0.000028;
	double twr5 = 0.000030;
	double twr6 = 0.000032;
	double twr7 = 0.000034;
	double twr8 = 0.000038;

	// for output
	double tprint = 0;
	ofstream outputFile("output.txt", ios::out);



/*****************************************************************************
 *
 * The subroutine "cshift1R" performs a circular shift over a distance
 * of 1 to the right on an array doubles of size n
 */ 
void cshift1R(double array[], const int n)
{	// it is supposed that the array has n elements
	// and it shifts over a distance of 1 "to the right"
	double last = array[n-1];
	for(int i=n-1; i>0; --i) {
		array[i] = array[i-1];
	}
	array[0] = last;
}

/*****************************************************************************
 *
 * The subroutine "cshift1L" performs a circular shift over a distance
 * of 1 to the left on an array of doubles of size n
 */ 
void cshift1L(double array[], const int n)
{	// it is supposed that the array has n elements
	// and it shifts over a distance of 1 "to the left"
	double first = array[0];
	for(int i=0; i<n-1; ++i) {
		array[i] = array[i+1];
	}
	array[n] = first;
}

/*****************************************************************************
 *
 * The subroutine "cshift" performs a circular shift over a distance
 * of d on an array of doubles of size n
 */ 
void cshift(double array[], int d, const int n)
{	// it is supposed that the array has n elements
	if (d<0) {
		for(int i=0; i>d; --i) {cshift1L(array,n);}
	}
	else {
		if (d>0) {
			for(int i=0; i<d; ++i) {cshift1R(array,n);}
		}
	}
}

/*****************************************************************************
 *
 * The subroutine "anyOverlap" checks the expression below
 */
bool anyOverlap(double r[],double hh[], double ht[], const int n)
{  
	// I am not sure wether only r[i] has to be checked with hh[i]
	// perhaps all possible r[i] and hh[j] have to be considered
	// in Fortran: any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht))
	for(int i=0; i<n; ++i) {
		if (fabs(r[i])<=(2.0*hh[i])) return true;
		if (fabs(r[i])<=(2.0*ht[i])) return true;
	}
	return false;
}


/*****************************************************************************
 *
 * The subroutine "gradww" computes the gradient of the interpolating kernal
 */
void gradww(double r[],double hh[],double ht[],int n,double gradw[])
{
	double q1 = 1.0;
	double q2 = 2.0;

	for(int i=0; i<n; ++i) { 
		q[i] = (fabs(r[i])/hh[i]);
	}

	for(i=0; i<n; ++i) { 
		if (q[i]<=q1) {
			gradw[i] = (1.5*(r[i]*r[i])/pow(hh[i],3));
			gradw[i] = gradw[i]-(2.0*fabs(r[i])/(hh[i]*hh[i]));
			gradw[i] = gradw[i]/hh[i];
		}
	}

	for(i=0; i<n; ++i) { 
		if ((q[i]>q1)&&(q[i]<=q2)) {
			gradw[i] = pow((2.0-q[i]),2);
			gradw[i] = -3.0*gradw[i]/hh[i];
			gradw[i] = gradw[i]/(6.0*hh[i]);
		}
	}

	for(i=0; i<n; ++i) { 
		if (r[i]!=0.0) {
			gradw[i] = gradw[i]*r[i]/fabs(r[i]);
		}
		else {
			gradw[i] = 0.0;
		}
	}

	for(i=0; i<n; ++i) { 
		if (q[i]>q2) { 
			gradw[i] = 0.0;
		}
		q[i] = fabs(r[i])/ht[i];
	}

	for(i=0; i<n; ++i) { 
		if (q[i]<q1) { 
			gradw2[i] = (1.5*(r[i]*r[i])/pow(ht[i],3));
			gradw2[i] = gradw2[i]-(2.0*fabs(r[i])/(ht[i]*ht[i]));
			gradw2[i] = gradw2[i]/ht[i];
		}
	}

	for(i=0; i<n; ++i) { 
		if ((q[i]>q1)&&(q[i]<=q2)) {
			gradw2[i] = pow((2.0-q[i]),2);
			gradw2[i] = -3.0*gradw2[i]/ht[i];
			gradw2[i] = gradw2[i]/(6.0*ht[i]);
		}
	}

	for(i=0; i<n; ++i) { 
		if (r[i]!=0.0) {
			gradw2[i] = gradw2[i]*r[i]/fabs(r[i]);
		}
		else {
			gradw2[i] = 0.0;
		}
	}

	for(i=0; i<n; ++i) { 
		if (q[i]>q2) {
			gradw2[i] = 0.0;
		}
		gradw[i] = (gradw[i]+gradw2[i])/2.0;
	}
}

/*****************************************************************************
 *
 * the subroutine "accelerations" computes the accelerations on each particle
 * as well as the rate of change of the specific energy (or entropy), u.
 */
void accelerations (double xx[],double vx[],double ax[],double au[],
					double hh[],double p[],double u[],double mass[],
					double den[],double divv[],double sound[],
					double maxmu[],int n,bool inorm,double alpha,
					double beta,bool ibulk,double alphab,
					double betab,double gamma,bool ientro,
					bool ientro2)
{
	for(int i=0; i<n; ++i) {
		// initialize the accelerations and du/dt to zero
		ax[i] = 0.0;     
		au[i] = 0.0;

		// store copies of parameters for circular shift
		xt[i]		= xx[i];
		vxt[i]		= vx[i];
		ht[i]		= hh[i];
		axt[i]		= ax[i];
		aut[i]		= au[i];
		pt[i]		= p[i];
		ut[i]		= u[i];
		masst[i]	= mass[i];
		dent[i]		= den[i];
		divvt[i]	= divv[i];
		soundt[i]	= sound[i];
		maxmu[i]	= 0.0;
		maxmut[i]	= 0.0;
		r[i]		= 0.0;
	}
	int ncount = 0;    // a counter

	// compute forces while any two particles are overlapping
	do { //while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		// shift the data one location

		cshift(xt,+1,n);
		cshift(vxt,+1,n);
		cshift(axt,+1,n);
		cshift(ht,+1,n);
		cshift(pt,+1,n);
		cshift(ut,+1,n);
		cshift(aut,+1,n);
		cshift(masst,+1,n);
		cshift(dent,+1,n);
		cshift(divvt,+1,n);
		cshift(soundt,+1,n);
		cshift(maxmut,+1,n);

		for(int i=0; i<n; ++i) {
			r[i] = xx[i]-xt[i];        // the distance between the particles
			h[i] = (ht[i]+hh[i])/2.0;  // their average smoothing length
			v[i] = vx[i]-vxt[i];       // their relative velocity
		
			// molecular viscosity
			if (inorm==1) {
				eta[i] = 0.1*(h[i]*h[i]);
				mu[i]  = h[i]*v[i]*r[i]/((r[i]*r[i])+eta[i]);
				if ((v[i]*r[i])>0.0) mu[i] = 0.0;
				if (fabs(mu[i])>maxmu[i]) maxmu[i] = fabs(mu[i]);
				if (fabs(mu[i])>maxmut[i]) maxmut[i] = fabs(mu[i]);
				cij[i] = (sound[i]+soundt[i])/2.0;
				denij[i] = (den[i]+dent[i])/2.0;
				artv[i] = (-alpha*mu[i]*cij[i])+(beta*(mu[i]*mu[i]));
				artv[i] = artv[i]/denij[i];
			}

			// bulk viscosity
			if (ibulk==1) {
				if (divv[i]>=0.0) {
					qi[i] = 0.0;
				}
				else {
					qi[i] = (alphab*hh[i]*den[i]*sound[i]*fabs(divv[i]))+
							(betab*(hh[i]*hh[i])*den[i]*(divv[i]*divv[i]));
				}
				if (divvt[i]>=0.0) {
					qj[i] = 0.0;
				}
				else {
					qj[i] = (alphab*ht[i]*dent[i]*soundt[i]*fabs(divvt[i]))+
				            (betab*(ht[i]*ht[i])*dent[i]*(divvt[i]*divvt[i]));
				}
				artv[i] = (qi[i]/(den[i]*den[i]))+(qj[i]/(dent[i]*dent[i]));
				if ((r[i]*v[i])>0.0) artv[i] = 0.0;
			}
		}

		// next compute the gradient of w, gradw
		gradww(r,hh,ht,n,gradw);

		for(i=0; i<n; ++i) {
			// the forces
			if (ientro2==0) {
				// the two symmetric forms 
				//pforce[i] = (p[i]/(den[i]*den[i]))+(pt[i]/(dent[i]*dent[i])); 
				// for the pressure force
				pforce[i] = (2.0*sqrt(p[i]*pt[i])/(den[i]*dent[i]));
				// add the art. visc. term and multiply by the gradient of the int.kernal
				pforce[i] = (pforce[i]+artv[i])*gradw[i]; 
				
				ax[i] = ax[i]-(masst[i]*pforce[i]);		// add to the running 
				axt[i] = axt[i]+(mass[i]*pforce[i]);    // totals
			}
			else {
				// the following is for the case when the entropic 
				// formulation of the force is used
				pforce[i] = ut[i]+((gamma-1.0)*u[i]);
				pforce[i] = pforce[i]*pow(den[i],(gamma-2.0));

				ax[i] = ax[i]-(masst[i]*(pforce[i]+artv[i])*gradw[i]);

				pforce[i] = u[i]+((gamma-1.0)*ut[i]);
				pforce[i] = pforce[i]*pow(dent[i],(gamma-2.0));

				axt[i] = axt[i]+(mass[i]*(pforce[i]+artv[i])*gradw[i]);
			}

			// du/dt
			if (ientro==0) {
				au[i] = au[i]+(masst[i]*pforce[i]*v[i]);
				aut[i] = aut[i]+(mass[i]*pforce[i]*v[i]);
			}
			 else {
				// the rate of change of the entropic function
				au[i] = au[i]+(masst[i]*artv[i]*v[i]*gradw[i]);
				aut[i] = aut[i]+(mass[i]*artv[i]*v[i]*gradw[i]);
			}
		}
		ncount++;

	} while (anyOverlap(r,hh,ht,n));

	// shift the temporary totals back to their starting points and add to get
	// the complete sum
	cshift(axt,-ncount,n);
	cshift(aut,-ncount,n);
	
	for(i=0; i<n; ++i) {
		ax[i] = ax[i]+axt[i];
		au[i] = 0.5*(au[i]+aut[i]);
	}
	cshift(maxmut,-ncount,n);

	for(i=0; i<n; ++i) {
		if (maxmut[i]>maxmu[i]) maxmu[i] = maxmut[i];
	}
}

/*****************************************************************************
 *
 * The subroutine "therm_diffusion" computes an additional thermal diffusion
 * term which can be added to the energy equation
 */
void therm_diffusion(double xx[],double vx[],double hh[],
					 double mass[],double den[],double p[],
					 double u[],double divv[],double sound[],
					 int n,double g1,double g2,double maxmu[],
					 double gamma,double thermdiff[])
{
	for(int i=0; i<n; ++i) {
		// thermal diffusion
		thermdiff[i] = 0.0;
		thermdifft[i] = thermdiff[i];
		xt[i] = xx[i];
		vxt[i] = vx[i];
		ht[i] = hh[i];
		pt[i] = p[i];
		ut[i] = u[i];
		masst[i] = mass[i];
		dent[i] = den[i];
		divvt[i] = divv[i];
		sound[i] = sqrt(gamma*p[i]/den[i]);
		soundt[i] = sound[i];
		maxmu[i] = 0.0;
		r[i] = 0.0;
	}
	int ncount = 0;
	
	do { //while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		cshift(xt,+1,n);
		cshift(vxt,+1,n);
		cshift(ht,+1,n);
		cshift(pt,+1,n);
		cshift(ut,+1,n);
		cshift(masst,+1,n);
		cshift(dent,+1,n);
		cshift(divvt,+1,n);
		cshift(soundt,+1,n);
		cshift(thermdifft,+1,n);
		
		for(i=0; i<n; ++i) {
			r[i] = xx[i]-xt[i];
			uu[i] = u[i]-ut[i];

			qqi[i] = (g1*hh[i]*den[i]*sound[i])+
					 (g2*den[i]*(hh[i]*hh[i])*(fabs(divv[i])-divv[i]));

			qqj[i] = (g1*ht[i]*dent[i]*soundt[i])+
					 (g2*dent[i]*(ht[i]*ht[i])*(fabs(divvt[i])-divvt[i]));

			qqi[i] = qqi[i]/den[i];
			qqj[i] = qqj[i]/dent[i];
			qqj[i] = (qqi[i]+qqj[i])/2.0;

			denij[i] = (den[i]+dent[i])/2.0;
		}

		// next compute the gradient of w, gradw
		gradww(r,hh,ht,n,gradw);

		for(i=0; i<n; ++i) {
			h[i] = (hh[i]+ht[i])/2.0;
			eta[i] = 0.1*(h[i]*h[i]);

			if (r[i]!=0.0) {
				thermdiff[i] = thermdiff[i]+(masst[i]*qqj[i]*uu[i]*r[i]
							   *gradw[i]/(denij[i]*((r[i]*r[i])+eta[i])));
				thermdifft[i] = thermdifft[i]-(mass[i]*qqj[i]*uu[i]*r[i]
								*gradw[i]/(denij[i]*((r[i]*r[i])+eta[i])));
			}
		}
		ncount++;
	}
	while (anyOverlap(r,hh,ht,n));

	cshift(thermdifft,-ncount,n);

	for(i=0; i<n; ++i) {
		thermdiff[i] = 2.0*(thermdiff[i]+thermdifft[i]);
	}
}
        
/*****************************************************************************
 *
 * The subroutine "push1" advances the particles' positions, velocities, and 
 * specific energies forward in time by half a time step
 */
void push1(double xx[],double u[],double vx[],double ax[],
		   double au[],double c1,double c2,int n)
{        
	for(int i=0; i<n; ++i) { 
		vx[i] = vx[i]+(c1*ax[i]);
		xx[i] = xx[i]+(c2*ax[i])+(c1*vx[i]);
		u[i]  = u[i]+(c1*au[i]);
	}
}
/*****************************************************************************
 *
 * The subroutine "benz_h1" advances the particles' smoothing length forward
 * in time by half a time step according to the algorithm of Benz
 */
void benz_h1(double hh[],double divv[],double c1,int n)
{
	for(int i=0; i<n; ++i) {
           hh[i] = hh[i]+(c1*hh[i]*divv[i]);
	}
}
/*****************************************************************************
 *
 * The subroutine "push2" is called after the accelerations have been computed
 * and advances the particles through a complete time step
 */
void push2(double xx[],double u[],double vx[],double ax[],
		   double au[],double c1,double c2,double delt,
		   double xo[],double vxo[],double uo[],double t,int n)
{

	if (t==0.0) {
		for(int i=0; i<n; ++i) {
			xx[i] = xx[i]+(c1*vx[i]);
			vx[i] = vx[i]+(c1*ax[i]);
			u[i]  = u[i]+(c1*au[i]);
		}  
	}
	else {
		for(int i=0; i<n; ++i) {
			vx[i] = vxo[i]+(delt*ax[i]);
			xx[i] = xo[i]+(delt*(vxo[i]+(c1*ax[i])));
			u[i]  = uo[i]+(delt*au[i]);
		}
	}
}

/*****************************************************************************
 *
 * The subroutine "benz_h2" advances the smoothing lengths through a complete
 * time step and is called in consort with "push2"
 */
void benz_h2(double hh[],double delt,double divv[],double ho[],double c1,int n, double t)
{

	if (t==0.0) {
		for(int i=0; i<n; ++i) {
			hh[i] = hh[i]+(c1*divv[i]*hh[i]);
		}
	}
	else {
		for(int i=0; i<n; ++i) {
			hh[i]=ho[i]+(delt*divv[i]*hh[i]);
		}
	}
}

/*****************************************************************************
 * 
 * The subroutine "ww" computes the 1-d value of the interpolating kernal
 * given a distance (r) between the particles and their smoothing lengths
 */
void ww(double r[],double hh[],double ht[],const int n,double w[])
{
	double q1	= 1.0;
	double q2	= 2.0;

	for(int i=0; i<n; ++i) { 

		q[i] = fabs(r[i])/hh[i];

		if (q[i]<=q1) {
			w[i] = (2.0/3.0)-(q[i]*q[i])+(0.5*pow(q[i],3));
			w[i] = w[i]/hh[i];
		}
		else {
			if ((q[i]>q1)&&(q[i]<=q2)) {
				w[i] = pow((2.0-q[i]),3);
				w[i] = w[i]/(6.0*hh[i]);
			}
			else { 
				w[i] = 0.0;
			}
		}

		q[i] = fabs(r[i])/ht[i];

		if (q[i]<=q1) {
			w2[i] = (2.0/3.0)-(q[i]*q[i])+(0.5*pow(q[i],3));
			w2[i] = w2[i]/ht[i];
		}
		else {
			if ((q[i]>q1)&&(q[i]<=q2)) {
				w2[i] = pow((2.0-q[i]),3);
				w2[i] = w2[i]/(6.0*ht[i]);
			}
			else {
				w2[i]=0.0;
			}
		}
		
		w[i] = (w[i]+w2[i])/2.0;

	}
}


/*****************************************************************************
 *
 * The subroutine "update_den" computes the density according the postions and
 * the smoothing lengths of the particles
 */  
void update_den(double xx[],double mass[],double den[],double hh[],
				int n)
{
	int		ncount	= 0;

	for(int i=0; i<n; ++i) {
		xt[i]		= xx[i];
		ht[i]		= hh[i];
		masst[i]	= mass[i];
		den[i]		= 0.0;
		dent[i]		= den[i];
		r[i]		= 0.0;
	}

	do { //Fortran: while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		for(int i=0; i<n; ++i) {// distance between particles
			r[i] = xx[i]-xt[i];
		}	

		ww(r,hh,ht,n,w);		// w
  
		for(i=0; i<n; ++i) {	// running total of the density
			den[i] += (masst[i]*w[i]);					
		}

		if (ncount>0) {				// running total of the
			for(int i=0; i<n; ++i) {// shifted density			
				dent[i] = dent[i]+(mass[i]*w[i]);
			}
		}
		
		ncount++;
  
		// shift the data
		cshift(xt,+1,n);
		cshift(ht,+1,n);
		cshift(masst,+1,n);
		cshift(dent,+1,n);

	} while (anyOverlap(r,hh,ht,n));

	cshift(dent,-ncount,n);
	for(i=0; i<n; ++i) { den[i] += dent[i];}
}

/*****************************************************************************
 *
 * The subroutine update_p computes the pressure for each particle assuming an
 * ideal gas
 */
void update_p(double hh[],double mass[],double u[],double den[],
			  double p[],double gamma,int n,bool ientro)
{
	for(int i=0; i<n; ++i) {
		if (ientro==0) {
			p[i] = (gamma-1.0)*u[i]*den[i];
		}
		else {
			p[i] = u[i]*pow(den[i],gamma);
		}
	}
}
/*****************************************************************************
 *
 * The subroutine new_delt computes a new time step
 */
void new_delt(double hh[],double divv[],double sound[],bool inorm,
			  double alpha,double beta,bool ibulk,double alphab,
			  double betab,double maxmu[],int n,double c1,double c2,
			  double cour,double delt)
{
	// in the case the molecular viscosity is used
	if ((inorm==1)||((inorm==0)&&(ibulk==0))) {
		for(int i=0; i<n; ++i) {
			delti[i] = (hh[i]*fabs(divv[i]))+sound[i]+
				(1.2*((alpha*sound[i])+(beta*maxmu[i])));
			delti[i] = hh[i]/delti[i];
		}
	}

	// in the case the bulk viscosity is used
	if (ibulk==1) {
		for(int i=0; i<n; ++i) {
			ht[i] = betab*hh[i]*fabs(divv[i]);
			if (divv[i]>=0.0) ht[i]=0.0;
			delti[i] = (hh[i]*fabs(divv[i]))+sound[i]+
				(1.2*((alphab*sound[i])+ht[i]));
			delti[i] = hh[i]/delti[i];
		}
	}

	delt = delti[0];
	for(int i=0; i<n; ++i) {if (delti[i]<delt) delt = delti[i];}
	
	delt	= cour*delt;
	c1		= delt/2.0;
	c2		= (delt*delt)/4.0;
}

/*****************************************************************************
 *
 * The subroutine "reflecting" set the boundry conditions.
 * This is done by using 50 guard particles on either end of the domain which 
 * are the mirro images of particles in the domain
 */
void reflecting(double xx[],double vx[],double xo[],double vxo[],
				double den[],double mass[],double p[],double u[],
				double uo[],double hh[],double ho[],int n)
{ 
	for(int i=0; i<50; ++i) {
		xx[i]	= -xx[100-i]+xx[50];
		vx[i]	= -vx[100-i];
		xo[i]	= -xo[100-i]+xo[50];
		vxo[i]	= -vxo[100-i];
		den[i]	= den[100-i];
		mass[i] = mass[100-i];
		p[i]	= p[100-i];
		u[i]	= u[100-i];
		uo[i]	= uo[100-i];
		hh[i]	= hh[100-i];
		ho[i]	= ho[100-i];
	}

	for(i=0; i<50; ++i) {
		xx[n-2-i]	= xx[n-52]+(xx[n-52]-xx[n-102+i]);
		vx[n-2-i]	= -vx[n-102+i];
		xo[n-2-i]	= xo[n-52]+(xo[n-52]-xo[n-102+i]);
		vxo[n-2-i]	= -vxo[n-102+i];
		den[n-2-i]	= -(-den[n-102+i]);
		mass[n-2-i]	= -(-mass[n-102+i]);
		p[n-2-i]	= -(-p[n-102+i]);
		u[n-2-i]	= -(-u[n-102+i]);
		hh[n-2-i]	= -(-hh[n-102+i]);
	}
}

/*****************************************************************************
 *
 * The subroutine "double_shock_tube" set up the initial conditions for the
 * the double shock tube problem of Woodward and Collela
 */
void double_shock_tube(double xx[],double vx[],double p[],double u[],
					   double den[],double mass[],double xmax,
					   double xmin,double hh[],int partno[],
					   double gamma,double hfac,int n,bool ientro)
{
    double h = hfac*(xmax-xmin)/(n-100);

	for(int i=0; i<n; ++i) {

		partno[i]	= i;
		hh[i]		= h;
		xx[i]		= (i-51)*(xmax-xmin)/(n-100);
		den[i]		= 1.0;
		mass[i]		= den[i]*(xmax-xmin)/(n-100);
		vx[i]		= 0.0;

		if (xx[i]<0.1) p[i]=1000.0;
		if ((xx[i]>=0.1)&&(xx[i]<0.9)) p[i]=0.01;
		if (xx[i]>=0.9) p[i]=100.0;

		u[i]		= p[i]/den[i];
		u[i]		= u[i]/(gamma-1.0);

		if (ientro==1) u[i]=p[i]/pow(den[i],gamma);
		// here density is mass per unit length
	}
}

/*****************************************************************************
 *
 * The subroutine "smooth" smooths a variable (here represented by u).
 */
void smooth(double xx[],double hh[],double mass[],double den[],
			const int n,double u[])
{

	int ncount = 0;

	for(int i=0; i<n; ++i) { 
		xt[i]		= xx[i];
		ht[i]		= hh[i];
		masst[i]	= mass[i];
		old_u[i]	= u[i];
		old_ut[i]	= u[i];
		u[i]		= 0.0;
		ut[i]		= 0.0;
		r[i]		= 0.0;
	}

	do { //Fortran: while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		for(int i=0; i<n; ++i) { r[i] = xx[i]-xt[i];}

		ww(r,hh,ht,n,w);

		for(i=0; i<n; ++i) { 
			u[i] = u[i]+(masst[i]*old_ut[i]*w[i]);
		}

		if (ncount>0) { ut[i] = ut[i]+(mass[i]*old_u[i]*w[i]);}

		ncount++;

		cshift(xt,+1,n);
		cshift(ht,+1,n);
		cshift(masst,+1,n);
		cshift(old_ut,+1,n);
		cshift(ut,+1,n);

	} while (anyOverlap(r,hh,ht,n));


	cshift(ut,-ncount,n);

	for(i=0; i<n; ++i) { 
		u[i] += ut[i];
		u[i] /= den[i];
	}
}

/*****************************************************************************
 *
 * The subroutine "div_v" computes the divergence of the velocity field
 */
void div_v(double xx[],double vx[],double hh[],double mass[],
		   double den[],int n,double divv[])
{
	for(int i=0; i<n; ++i) {
		divv[i]		= 0.0;
		xt[i]		= xx[i];
		vxt[i]		= vx[i];
		ht[i]		= hh[i];
		divvt[i]	= 0.0;
		masst[i]	= mass[i];
		r[i]		= 0.0;
	}
	int ncount = 0;

	do { //while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		cshift(xt,+1,n);
		cshift(vxt,+1,n);
		cshift(ht,+1,n);
		cshift(divvt,+1,n);
		cshift(masst,+1,n);

		for(int i=0; i<n; ++i) {
			r[i] = xx[i]-xt[i];
			v[i] = vx[i]-vxt[i];
		}
		// next compute the gradient of w, gradw
		gradww(r,hh,ht,n,gradw);

		for(i=0; i<n; ++i) {
			divv[i] = divv[i]+(masst[i]*v[i]*gradw[i]);
			// v is in the opposite sense as well as gradw, hence the +
			divvt[i] = divvt[i]+(mass[i]*v[i]*gradw[i]);
		}

		ncount++;

	} while (anyOverlap(r,hh,ht,n));

	cshift(divvt,-ncount,n);

	for(i=0; i<n; ++i) {
		divv[i] = -(divv[i]+divvt[i])/den[i];
	}
}

/*****************************************************************************
 *
 * The subroutine "smooth_den" computes a smoothed value for the density
 * for the Steinmetz and Muller alg.
 */
void smooth_den(double xx[],double hh[],double mass[],double den[],
				int n,double densmooth[])
{
	// compute a new smoothed density
	int ncount = 0;
	double hfac = 1.5;

	for(int i=0; i<n; ++i) {
		xt[i]			= xx[i];
		ht[i]			= hh[i];
		masst[i]		= mass[i];
		dent[i]			= den[i];
		densmooth[i]	= 0.0;
		denst[i]		= densmooth[i];
		r[i]			= 0.0;
	}

	do { //while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		for(int i=0; i<n; ++i) {
			r[i] = xx[i]-xt[i];
		}

		ww(r,hh,ht,n,w);

		for(i=0; i<n; ++i) {
			densmooth[i] += masst[i]*w[i]*dent[i];
		}

		if (ncount>0) {
			for(int i=0; i<n; ++i) {
				denst[i] += mass[i]*w[i]*den[i];
			}
		}

		ncount++;

		cshift(xt,+1,n);
		cshift(ht,+1,n);
		cshift(masst,+1,n);
		cshift(dent,+1,n);
		cshift(denst,+1,n);

	} while (anyOverlap(r,hh,ht,n));

	cshift(denst,-ncount,n);
	for(i=0; i<n; ++i) {
		densmooth[i] = (densmooth[i]+denst[i])/den[i];
		hh[i] = hfac*(mass[i]/densmooth[i]);
	}
}

/*****************************************************************************
 *
 * The subroutine "d_smooth_den" computes the rate of change of the smoothed
 * density for the Steinmetz + Muller alg.
 */
void d_smooth_den(double xx[],double hh[],double mass[],double den[],
				  double densmooth[],double divv[],int n,
				  double dsmoothdt[])
{
	int ncount = 0;

	for(int i=0; i<n; ++i) {
		xt[i] = xx[i];
		ht[i] = hh[i];
		masst[i] = mass[i];
		dent[i] = densmooth[i];
		divvt[i] = divv[i];
		dsmoothdt[i] = 0.0;
		dsmt[i] = 0.0;
		r[i] = 0.0;
	}

	do { //while (any(abs(r).le.(2.0*hh).or.abs(r).le.(2.0*ht)))

		for(int i=0; i<n; ++i) {
			r[i] = xx[i] - xt[i];
		}

		ww(r,hh,ht,n,w);

		for(i=0; i<n; ++i) {
			dsmoothdt[i] += masst[i]*dent[i]*divvt[i]*w[i];
		}
		
		if (ncount>0) {
			dsmt[i] += mass[i]*densmooth[i]*divv[i]*w[i];
		}
		ncount++;

		cshift(xt,+1,n);
		cshift(ht,+1,n);
		cshift(masst,+1,n);
		cshift(dent,+1,n);
		cshift(divvt,+1,n);
		cshift(dsmt,+1,n);

	} while (anyOverlap(r,hh,ht,n));

	cshift(dsmt,-ncount,n);
	for(i=0; i<n; ++i) {
		dsmoothdt[i] += dsmt[i];
		dsmoothdt[i] = 2.0*dsmoothdt[i]/den[i];
		dsmoothdt[i] = (densmooth[i]*divv[i])-dsmoothdt[i];
	}
}
       
/*****************************************************************************
 *
 * The subroutine SM_h computes new smoothing lengths according to the alg. of
 * Steinmetx and Muller
 */
void SM_h(double xx[],double hh[],double mass[],double den[],
		  double divv[],double t,double delt,int n)
{
	smooth_den(xx,hh,mass,den,n,densmooth); // compute a smoothed 
											// value of the density

	d_smooth_den(xx,hh,mass,den,densmooth,divv,n,dsmoothdt);
											// compute the rate of
											// of change of the 
											// smoothed density
	// 1, estimate the smoothed density
	if (t==0.0) for(int i=0; i<n; ++i) {densmooth[i] = den[i];}
	for(int i=0; i<n; ++i) {densmooth[i] += (delt*dsmoothdt[i]);}

	// 2, calculate the average density of the system
	double sum = 0;
	for(i=0; i<n; ++i) {
		sum += 1.0/densmooth[i];
	}
	double avden = 1.0/sum;

	// 3, find new smoothing length
	double hfac = 1.5;
	for(i=0; i<n; ++i) {avh[i] = hfac*mass[i]/avden;}

	// for 3-d need to use (mass/den)**1/3 and (mass/den)**1/2 for 2-d
	int kkk = 1;
	for(i=0; i<n; ++i) {
		hh[i] = avh[i]*pow((avden/densmooth[i]),kkk);
	}
}
//*****************************************************************************


void main() {

	/****************************************************************************
	 *
	 *
	 *^^^^^^^^^^^^^^^^^^^^^^BEGIN INITIALIZATIONS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 */	
	t		= 0.0;				// start with time zero
	delt	= 1.0e-6;			// the initial value of the time step is set 
								// to a small value, it is modified as the 
								// calculation procedes

    c1		= delt/2.0;			// c1 and c2 are constants used in the integration of
    c2		= (delt*delt)/4.0;	// the equation of motion 

    tend	= 0.00003801;			// the ending time is set

    gamma	= 1.4;				// the ratio of specific heats is set

    cour	= 0.3;				// the courant number is set

	nsteps = 0;					// no steps taken when we start a new simulation

	xmin = 0.0;					// the minimum and maximum values of x, the define the 
	xmax = 1.0;					// domain of the problem

    alpha=1.0;					// alpha and beta are the coefficients used in the 
    beta=2.0;					// molecular viscosity

    alphab=0.5;					// alphab and betab are the coefficients used in the 
    betab=1.5;					// bulk viscosity

    g1=1.0/10.0;				// g1 and g2 are the coefficients used in the 
    g2=1.0;						// thermal diffusion term

    hfac=2.0;					// hfac is a factor which determines the initial 
								// smoothing length of particles, 
								// hh = hfac * pow((mass/den),(1/d)), where d is the 
								// dimensionality
	tote=0.0;					// start with zero energy

	if (!outputFile){cout<<"Error can't open file 'output.txt' for output.\n";}

	/*^^^^^^^^^^^^^^^^^^^^END OF INITIALIZATIONS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 *
	 ******************************************************************************
	 *
	 *^^^^^^^^^^^^^^^^^^^SET UP OF INITIAL CONDITIONS^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 *
	 *     The following set of subroutines can be called individually to set up
	 *     the initial conditions for some different problems.
	 *
	 */
	
	// "double_shock_tube" sets up the initial conditions for the double shock tube
	// problem of Woodward and Collela
	double_shock_tube(xx,vx,p,u,den,mass,xmax,xmin,hh,partno,gamma,hfac,n,ientro);

	// "sod" sets up the initial conditions for Sod's problem
//    sod(xx,vx,mass,den,p,u,hh,gamma,xmax,xmin,partno,hfac,n,ientro);

	// "sod_const" is used if the distance between particles is constant
//	sod_const(xx,vx,mass,den,p,u,hh,gamma,xmax,xmin,partno,hfac,n,ientro);

	// "noh" sets up the initial conditions for the "shocks of noh", i.e. two 
	// colliding streams of gas, note: turn off reflecting boundries when doing this
	// problem
//	noh(xx,vx,mass,den,p,u,hh,gamma,xmax,xmin,partno,hfac,n,ientro);


	/******************************************************************************
	 *
	 *^^^^^^^^^^^^^^^^^^INITIAL SMOOTHING^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 *
	 *     Initial smoothing of the pressure (energy) and velocity fields is 
	 *     performed following Hernquist and Katz (1989).
	 */    
//Jim//	update_den(xx,mass,den,hh,n);

	//     smoothing the initial energy (entropy) field, u
//Jim//	smooth(xx,hh,mass,den,n,u);

	//     smoothing the initial velocity field, vx
//Jim//    smooth(xx,hh,mass,den,n,vx);

/*Jim	for(int i=0; i<n; ++i) {
		
		p[i]=(gamma-1.0)*u[i]*den[i];   // computes the pressure via the equation of state
		
		if (ientro==1) { p[i]=u[i]*pow(den[i],gamma);}
	}
*/	  
	if (ientro==1){
		
		for(int i=0; i<n; ++i) { p[i]=u[i];}
		
		smooth(xx,hh,mass,den,n,p);		// if the entropic formulation is used
										// compute the pressure correspondingly
		
		for(i=0; i<n; ++i) { p[i]=p[i]*pow(den[i],gamma);}
	}
	  
	/******************************************************************************
	 *
	 *^^^^^^^^^^^^^^^^^INITIAL SORT OF PARTICLES^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 *
	 *     The particles are sorted on x, this is not really necessary for one
	 *     dimensional problems since the particles do not become disordered in the
	 *     x coordinate.
	 */    
	r4_sort(xx,isort,partno,n);

	// Fortran: vx=vx(isort) etc.
	// sorting operation
	for(int i=0; i<n; ++i) {
		vxTemp[i]	= vx[i];
		axTemp[i]	= ax[i];
		pTemp[i]	= p[i];
		denTemp[i]	= den[i];
		hhTemp[i]	= hh[i];
		massTemp[i]	= mass[i];
		uTemp[i]	= u[i];
		auTemp[i]	= au[i];
	}

	for(i=0; i<n; ++i) {
		vx[i]	= vxTemp[isort[i]];
		ax[i]	= axTemp[isort[i]];
		p[i]	= pTemp[isort[i]];
		den[i]	= denTemp[isort[i]];
		hh[i]	= hhTemp[isort[i]];
		mass[i]	= massTemp[isort[i]];
		u[i]	= uTemp[isort[i]];
		au[i]	= auTemp[isort[i]];
	}


	/*^^^^^^^^^^^^^^^^^ALL INITIAL COMPUTATIONS NOW PERFORMED^^^^^^^^^^^^^^^^^^^^^^
	 *
	 *
	 ******************************************************************************
	 *
	 *
	 *^^^^^^^^^^^^^^^^^^^^THE MAIN TIME LOOP NOW BEGINS^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 */

	// setting the timer to zero
	TimeDelta();

    do { // do until the time (t) exceeds the ending time (tend)

		/************************************************************************** 
		 *
		 *     STORE the values of xx,vx,u and hh at the begining of the time step
		 *     to be used later in the particle pushing routine "push2"
		 */
		for(int i=0; i<n; ++i) {
			xo  [i]	= xx [i];
			vxo [i]	= vx [i];
			uo  [i]	= u  [i];
			ho  [i]	= hh [i];
			deno[i]	= den[i];
		}

		/**************************************************************************
		 *
		 *     PUSHING THE PARTICLES 1/2 A TIME STEP AHEAD OF THE ACCELARATIONS
		 */
		SM_h(xx,hh,mass,den,divv,t,delt,n); // update values of hh
                                            // according to 
                                            // Steinmetz and Muller


		if (t==0.0) 
			goto sortIt;     // at the first time step, 
											// skip the 1/2 step

		push1(xx,u,vx,ax,au,c1,c2,n);       // push the particles 1/2 a time
											// step ahead

//      benz_h1(hh,divv,c1,n);				// advance hh 1/2 a time step ahead
											// according to Benz's alg.

		/**************************************************************************
		 *
		 *     SORTING
		 */
sortIt:
		r4_sort(xx,isort,partno,n); //sort on x after push

		// Fortran: vx=vx(isort) etc.
		// sorting operation
		for(i=0; i<n; ++i) {
			vxTemp[i]	= vx[i];
			axTemp[i]	= ax[i];
			pTemp[i]	= p[i];
			denTemp[i]	= den[i];
			hhTemp[i]	= hh[i];
			massTemp[i]	= mass[i];
			uTemp[i]	= u[i];
			auTemp[i]	= au[i];
			xoTemp[i]	= xo[i];
			vxoTemp[i]	= vxo[i];
			uoTemp[i]	= uo[i];
			hoTemp[i]	= ho[i];
		}

		for(i=0; i<n; ++i) {
			vx[i]	= vxTemp[isort[i]];
			ax[i]	= axTemp[isort[i]];
			p[i]	= pTemp[isort[i]];
			den[i]	= denTemp[isort[i]];
			hh[i]	= hhTemp[isort[i]];
			mass[i]	= massTemp[isort[i]];
			u[i]	= uTemp[isort[i]];
			au[i]	= auTemp[isort[i]];
			xo[i]	= xoTemp[isort[i]];
			vxo[i]	= vxoTemp[isort[i]];
			uo[i]	= uoTemp[isort[i]];
			ho[i]	= hoTemp[isort[i]];
		}

		/**************************************************************************
		 *
		 *     SETTING BOUNDRY CONDITIONS, subroutine "reflecting" creates
		 *     reflecting (solid wall) boundry conditions
		 */
		if (irefl==1) {
			reflecting(xx,vx,xo,vxo,den,mass,p,u,uo,hh,ho,n);
		}

		/**************************************************************************
		 *
		 * the total energy
		 */
/*Jim//		if (t==0.0) {
			for(int i=0; i<n; ++i) {
				vxt[i]=vx[i];
				vxt[i]=0.5*vxt[i]*vxt[i];
			}
			if (ientro==0) {
				if (irefl==1) {
					//toteo=sum(mass*u,mask=partno.gt.51.and.partno.lt.(n-50));
					double sum = 0;
					for(i=52; i<(n-50); ++i) {
						sum += mass[i]*u[i];
					}
					toteo=sum;
				}
				else {
					//toteo=sum(mass*u);
					double sum = 0;
					for(i=0; i<n; ++i) {
						sum += mass[i]*u[i];
					}
					toteo=sum;
				}
			}
			else {
				if (irefl==1) {
					//toteo=sum(mass*u*(den**(gamma-1.0))/(gamma-1.0),mask=partno.gt.51.and.partno.lt.(n-50));
				    double sum = 0;
				    for(i=52; i<(n-50); ++i) {
						sum += mass[i]*u[i]*(pow(den[i],(gamma-1.0))/(gamma-1.0));
					}
					toteo=sum;
				}
				else {
					//toteo=sum(mass*u*(den**(gamma-1.0))/(gamma-1.0));
					double sum = 0;
					for(i=0; i<n; ++i) {
						sum += mass[i]*u[i]*(pow(den[i],(gamma-1.0))/(gamma-1.0));
					}
					toteo=sum;
				}
			}

			if (irefl==1) {
				//toteo=toteo+sum(mass*vxt,mask=partno.gt.51.and.partno.lt.(n-50));
				double sum = 0;
				for(i=52; i<(n-50); ++i) {
					sum += mass[i]*vxt[i];
				}
				toteo+=sum;
			}
			else {
				//toteo=(toteo+sum(mass*vxt));
			   	double sum = 0;
				for(i=0; i<n; ++i) {
					sum += mass[i]*vxt[i];
				}
				toteo+=sum;
			}
		}
//Jim*/
		/**************************************************************************
		 *
		 *     SOUND SPEED 
		 */
		for(i=0; i<n; ++i) { sound[i]=sqrt(gamma*p[i]/den[i]);}

		/**************************************************************************
		 *
		 *     THERMAL DIFFUSION
		 *
		 *     Here, the subroutine "therm_diffusion" is used, which computes
		 *     an additional term to be added to the energy equation for each
		 *     particle as described by Monaghan (1988), unpublished preprint.
		 */
		if (itherm==1) {    // if the switch is set, compute it
			therm_diffusion(xx,vx,hh,mass,den,p,u,divv,sound,n,g1,g2,maxmu,gamma,thermdiff); 
		}
      
		/******************************************************************************
		 *
		 *     ACCELERATIONS AND DU/DT
		 *
		 *     The subroutine "accelerations" computes the acceleration of each
		 *     particle (ax) due to pressure forces and artificial viscosity as well as
		 *     heating of the particles (au = du/dt).
		 */
		accelerations(xx,vx,ax,au,hh,p,u,mass,den,divv,
			sound,maxmu,n,inorm,alpha,beta,ibulk,alphab,betab,
			gamma,ientro,ientro2);
      
		if (ientro==1) { 
			for(i=0; i<n; ++i) { 
				au[i]=au[i]*u[i]*den[i]*(gamma-1.0)/p[i];
			}
		}

		//**************************************************************************

		if (itherm==1) { // if the switch is set, add the thermal
                         // diffusion term
			for(i=0; i<n; ++i) { au[i]=au[i]+thermdiff[i];}              
         
		}

		/**************************************************************************
		 *        
		 *     PUSH THE PARTICLES the rest of the way through the time step 
		 */
		push2(xx,u,vx,ax,au,c1,c2,delt,xo,vxo,uo,t,n);


		/**************************************************************************
		 *
		 *     UPDATE HH, following Benz (1990)
		 */        
		//      benz_h2(hh,delt,divv,ho,c1,n,t);
        
		/**************************************************************************
		 *
		 *     SORTING into x order after pushing particles
		 */
		r4_sort(xx,isort,partno,n);
      
		for(i=0; i<n; ++i) { 
			vxTemp[i]	= vx[i];
			axTemp[i]	= ax[i];
			pTemp[i]	= p[i];
			denTemp[i]	= den[i];
			hhTemp[i]	= hh[i];
			massTemp[i]	= mass[i];
			uTemp[i]	= u[i];
			auTemp[i]	= au[i];
		}
        
		for(i=0; i<n; ++i) { 
			vx[i]	= vxTemp[isort[i]];
			ax[i]	= axTemp[isort[i]];
			p[i]	= pTemp[isort[i]];
			den[i]	= denTemp[isort[i]];
			hh[i]	= hhTemp[isort[i]];
			mass[i] = massTemp[isort[i]];
			u[i]	= uTemp[isort[i]];
			au[i]	= auTemp[isort[i]];
		}

		/**************************************************************************
		 *      
		 *     UPDATE THE DENSITY AND PRESSURE with the new positions and new values
		 *     of the energy density.
		 */
		update_den(xx,mass,den,hh,n);

		if (ientro==1) {

			for(i=0; i<n; ++i) { p[i]=u[i];}
			smooth(xx,hh,mass,den,n,p);  // if the entropic formulation is used
                                         // compute the pressure correspondingly
			
			for(i=0; i<n; ++i) { p[i]=p[i]*pow(den[i],gamma);}
		}
		else {

			update_p(hh,mass,u,den,p,gamma,n,ientro);
		}

		//where (p.le.0.0) p=0.0
		for(i=0; i<n; ++i) {
			if(p[i]<=0.0) { p[i]=0.0;}
		}

		for(i=0; i<n; ++i) { 
			sound[i]=sqrt(gamma*p[i]/den[i]);	// a new sound speed is 
		}										// computed for each particle

		/**************************************************************************
		 *
		 * total energy
		 */        
/*Jim//		for(i=0; i<n; ++i) { 
			vxt[i]=vx[i];
			vxt[i]=0.5*vxt[i]*vxt[i];
		}

        if (ientro==0) {
			if (irefl==1) {
				//tote=sum(mass*u,mask=partno.gt.51.and.partno.lt.(n-50))
				double sum = 0;
				for(i=52; i<(n-50); ++i) {
					sum += mass[i]*u[i];
				}
				toteo=sum;
			}
			else {
				//tote=sum(mass*u)
				double sum = 0;
				for(i=0; i<n; ++i) {
					sum += mass[i]*u[i];
				}
				toteo=sum;
			}
		}
        else {
			if (irefl==1) {
				//tote=sum(mass*u*(den**(gamma-1.0))/(gamma-1.0),mask=partno.gt.51.and.partno.lt.(n-50))
				double sum = 0;
				for(i=52; i<(n-50); ++i) {
					sum += mass[i]*u[i]*(pow(den[i],(gamma-1.0))/(gamma-1.0));
				}
				toteo=sum;
			}
			else {
				//tote=sum(mass*u*(den**(gamma-1.0))/(gamma-1.0))
				double sum = 0;
				for(i=0; i<n; ++i) {
					sum += mass[i]*u[i]*(pow(den[i],(gamma-1.0))/(gamma-1.0));
				}
				toteo=sum;
			}
        }

        if (irefl==1) {
			//tote=tote+sum(mass*vxt,mask=i.gt.51.and.partno.lt.(n-50))
			double sum = 0;
			for(i=52; i<(n-50); ++i) {
				sum += mass[i]*vxt[i];
			}
			toteo+=sum;
		}
		else {
			//toteo=(toteo+sum(mass*vxt));
		  	double sum = 0;
			for(i=0; i<n; ++i) {
			    sum += mass[i]*vxt[i];
			}
			toteo+=sum;
		}
		
		double  perc=fabs(toteo-tote)/toteo;
        perc=100.0*perc;

//      print *,' tote = ',tote,' perc = ',perc
//Jim*/
		/**************************************************************************
		 *
		 *     BOUNDRY conditions are reset since the paticles have been pushed a 
		 *     second time
		 */        
		if (irefl==1) {
			reflecting(xx,vx,xo,vxo,den,mass,p,u,uo,hh,ho,n);
		}

		/**************************************************************************
		 *
		 *     computing the DIVERGENCE OF THE VELOCITY FIELD following Monaghan (1988)
		 *     unpublished preprint.
		 */      
		div_v(xx,vx,hh,mass,den,n,divv);
      
		/**************************************************************************
		 *
		 *     A NEW TIME STEP IS COMPUTED according to the courant condition following
		 *     Hernquist and Katz (1989), also Monaghan (1988).
		 */
		new_delt(hh,divv,sound,inorm,alpha,beta,ibulk,alphab,betab,maxmu,n,c1,c2,cour,delt);


		//**************************************************************************

		nsteps++;

		/**************************************************************************
		 *
		 *     OUTPUT
		 */
		bool timeToPrint = false;
		if ((t>twr1)&&(t<=(twr1+delt))) timeToPrint = true;
		if ((t>twr2)&&(t<=(twr2+delt))) timeToPrint = true;
		if ((t>twr3)&&(t<=(twr3+delt))) timeToPrint = true;
		if ((t>twr4)&&(t<=(twr4+delt))) timeToPrint = true;
		if ((t>twr5)&&(t<=(twr5+delt))) timeToPrint = true;
		if ((t>twr6)&&(t<=(twr6+delt))) timeToPrint = true;
		if ((t>twr7)&&(t<=(twr7+delt))) timeToPrint = true;
		if ((t>twr8)&&(t<=(twr8+delt))) timeToPrint = true;
		if (timeToPrint) {
			for(i=0; i<n; i++) {
				outputFile<<xx[i]<<" "<<den[i]<<" "<<p[i]<<" "<<vx[i]<<" "<<u[i]<<" "<<hh[i]<<endl;
			}
		}
        
		//*************************************************************************

		t=t+delt; // increment the time

		if((t>tprint)&&(t<=(tprint+delt))){
			cout<<"Time is: "<<t<<" completed in: "<<TimeDelta()<<" s"<<endl; cout.flush(); // let the user know something is happening
			tprint+=0.000001;
		}
	} while (t<=tend);

	/*^^^^^^^^^^^^^^^^^^^END OF TIME LOOP^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	 *
	 *
	 ******************************************************************************
	 */

} // end of main()

/*^^^^^^^^^^^^^^^^^END OF PROGRAM SPH-1D^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 *
 **********************************************************************************
 */