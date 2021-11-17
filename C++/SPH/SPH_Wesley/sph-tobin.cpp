/*********************************************************************
 * Translated from sph.f in Fortran to C++ by Wesley Miaw
 * <wesley@wesman.net>. Original code by Kevin Olson
 * <olson@jeans.gsfc.nasa.gov> is available from NASA/Goddard at
 * http://sdcd.gsfc.nasa.gov/ESS/exchange/contrib/olson/sph_1d.html.
 *
 * Included below is the original readme text from his sph.f code.
 *
 * The one dimensional code has been modified to work in three
 * dimensions.
 *
 * Modified: 05/08/2001
 *********************************************************************
 *
 *                             SPH_1D
 *
 *                           KEVIN OLSON
 *                       USRA and NCCS, GSFC
 *
 *                          April, 1993
 *
 * Smooth Particle Hydrodynamics is a gridless particle method for
 * solving the Lagrangian equations of hydrodynamics.  This
 * implementation includes variable smoothing lengths using both the
 * algorithm described in Benz (1990) and that described in Steinmetz
 * and Muller (1992).  Also included is a routine to compute an
 * additional thermal diffusion term which is added to the energy
 * equation (Monaghan 1988).  I have found that in problems where
 * the initial conditions are discontinuous, large deviations in the
 * energy profile occur at the contact discontinuity which forms near
 * the initial discontinuity.  Addition of a thermal diffusion term
 * removes this problem for the double shock tube problem (Woodward
 * and Colella, 1984).  Also included in the model are routines for
 * both the molecular and bulk viscosities.
 *
 * This implementation of SPH has been written in Maspar Fortran to
 * run on the Maspar MP-1 at GSFC.
 *
 * The following references where used extensively in building this
 * code:
 *
 *  Benz, W., 1990, in Numerical Modelling of Nonlinear Stellar
 *     Pulsations, Problems and Prospects, ed. J.R. Buchler, Kulwer
 *     Acdemics Publishers, Dordecht, Boston, London
 *
 *  Hernquist, L., and Katz, N., 1989, ApJ. Suppl., Vol. 70, 419
 *
 *  Monaghan, J.J., 1985, Comp. Phys. Rep., Vol. 3, 71
 *
 *  Monaghan, J.J., 1988, "SPH Meets the Shocks of Noh", unpublished
 *     preprint
 *
 *  Steinmetz, M., and Muller, E., 1992, preprint sub. to A&A
 *
 *  Woodward, P., and Colella, P., 1984, J. Comp. Phys., Vol. 54, 115
 *
 ********************************************************************/
#include "r4_sort.h"
#include "sph_mods.h"
#include "sod_int.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int debug = 2;


void main2(int argv, char** argc) {

// The number of particles, plus 100 for the boundaries.
static const int nsize = 400 + 100;

// Some global variables.
float xo[nsize], vxo[nsize], uo[nsize], vt[nsize], deno[nsize];
float yo[nsize], vyo[nsize];
float zo[nsize], vzo[nsize];
float alpha, beta, c1, c2, tend, alphab, betab;
float g1, g2, hfac, b1, grav;
float ho[nsize], sound[nsize];
float maxmu_x[nsize], maxmu_y[nsize], maxmu_z[nsize];
float tote, toteo;
float xs[128][128], ys[128][128];
int isort[nsize];

// Variable list, these definitions retain their meanings across
// subroutine boundaries.

// xx, vx, ax: position, velocity, and acceleration of a particle.
// yy, vy, ay: y-values.
// zz, vz, az: z-values.
float xx[nsize], vx[nsize], ax[nsize];
float yy[nsize], vy[nsize], ay[nsize];
float zz[nsize], vz[nsize], az[nsize];

// u, au: energy per unit mass and its rate of change, du/dt.
float u[nsize], au[nsize];

// den, p, mass: density, pressure, and mass of a particle.
float den[nsize], p[nsize], mass[nsize];

// hh: the smoothing length of a particle.
float hh[nsize];

// t, delt: the time and the time step, delta-t.
float t, delt;

// divv: the divergence of v.
float divv_x[nsize], divv_y[nsize], divv_z[nsize];

// thermdiff: the additional thermal diffusion term added to du/dt.
float thermdiff[nsize];

// gamma: the ratio of specific heats.
float gamma;

// cour: the courant number.
float cour;

// partno: the number of a particle in the list.
int partno[nsize];

/*********************************************************************
 * Begin Initializations
 ********************************************************************/

 if (debug>0) printf("begin initializations\n");

int n = nsize;       /* The total number of particles. */

delt = 1.0e-6;       /* The initial value of the time step is set to a
			small value, it is modified as the calculation
			procedes. */

c1 = delt / 2.0;     /* c1 and c2 are constants used in the */
c2 = pow(delt,2) / 4.0;  /*	integration of the equation of motion. */

tend = 0.035;  /* The ending time. */

gamma = 1.4;         /* The ratio of specific heats is set. */

cour = 0.3;          /* The courant number is set. */

int nsteps = 0;      /* A counter which counts the number of time
			steps taken. */

float xmin = 0.0;    /* The minimum and maximum values of x, y, and */
float xmax = 1.0;    /* z, they define the domain of the problem. */
float ymin = 0.0;
float ymax = 1.0;
float zmin = 0.0;
float zmax = 1.0;

alpha = 1.0;         /* alpha and beta are the coefficients used in */
beta = 2.0;          /* the molecular viscosity. */

alphab = 0.5;        /* alphab and betab are the coefficients used */
betab = 1.5;         /* in the bulk viscosity. */

g1 = 1.0 / 10.0;     /* g1 and g2 are the coefficients used in the */
g2 = 1.0;            /* thermal diffusion term. */


b1 = 0.2;            /* b1 is the coefficient used in the bouyancy
			calculation. It is multiplied by the
			difference in density. */

grav = 1.0;          /* grav is added to the bouyancy coefficient when
			cacluating the downward effect of
			bouyancy. This takes into account the effect
			of gravity on a denser particle falling
			down. */

hfac = 2.0;          /* hfac is a factor which determines the initial
			smoothing length of particles: hh = hfac *
			(mass / den)**(1/d), where d is the
			dimensionality. */

/*********************************************************************
 * The following switches are to be set to determine if and which
 * artificial viscosity is to be used, as well as if the thermal
 * diffusion term in the energy equation is to be included. A value of
 * 1 indicates that the switch is on.
 ********************************************************************/
int ibulk = 0;       /* The bulk viscosity switch. */

int inorm = 1;       /* The molecular viscosity switch. */

int itherm = 0;      /* The thermal diffusion switch. Note: do not use
			thermal diffusion with entropy. */

int ientro = 0;      /* Specifies whether or not the entropic equation
			is to be used in place of the energy
			equation. */

int ientro2 = 0;     /* Specifies whether entropy is used in lieu of
			the force equation. */

int irefl = 1;       /* Specifies whether reflecting boundaries are
			used. */

/*********************************************************************
 * Set up the initial conditions.
 ********************************************************************/

 if (debug > 0) printf("setting up initial conditions\n");
 
 sod(xx, vx, yy, vy, zz, vz, mass, den, p, u, hh, gamma, 
	 xmax, xmin, partno, hfac, n, ientro);
/*********************************************************************
 * Perform the initial smoothing.
 ********************************************************************/

 // Initial smoothing of the pressure (energy) and velocity fields is
 // performed following Hernquist and Katz (1989).

 if (debug > 0) printf("calling update_den\n");
 update_den(xx, yy, zz, mass, den, hh, n);

 // Smoothing the initial energy (entropy) field, u.
 
 if (debug > 0) printf("smooth(xx, yy, zz, hh, mass, den, n, u)\n");
 smooth(xx, yy, zz, hh, mass, den, n, u);

 // Smoothing the initial velocity field, vx.

 if (debug > 0) printf("smooth(xx, yy, zz, hh, mass, den, n, vx)\n");
 smooth(xx, yy, zz, hh, mass, den, n, vx);

 // Compute the pressure via the equation of state.
 if (debug > 0) printf("computing pressure via equation of state\n");
 for(int i = 0; i < n; i++) {
   p[i] = (gamma - 1.0) * u[i] * den[i];
   if (ientro == 1) p[i] = u[i] * pow(den[i],gamma);
 }

 if (ientro == 1) {
   for(int i = 0; i < n; i++) {
     p[i] = u[i];
   }
   // If the entropic formulation is used compute the pressure
   // accordingly.
   smooth(xx, yy, zz, hh, mass, den, n, p);
   for(i = 0; i < n; i++) {
     p[i] *= pow(den[i],gamma);
   }
 }
 
/*********************************************************************
 * Initial sort of particles.
 ********************************************************************/
// The particles are sorted on x, this is not really necessary for one
// dimensional problems since the particles do not become disordered
// in the x coordinate.
if (debug>0) printf("performing initial sort of particles\n");
r4_sort(xx, isort, partno, n);
float *swap_vx = new float[n], *swap_ax = new float[n];
float *swap_yy = new float[n], *swap_vy = new float[n], *swap_ay = new float[n];
float *swap_zz = new float[n], *swap_vz = new float[n], *swap_az = new float[n];
float *swap_p = new float[n], *swap_den = new float[n], *swap_hh = new float[n], *swap_mass = new float[n], *swap_u = new float[n], *swap_au = new float[n];
for(i = 0; i < n; i++) {
    swap_vx[i] = vx[isort[i]];
    swap_ax[i] = ax[isort[i]];
    swap_yy[i] = yy[isort[i]];
    swap_vy[i] = vy[isort[i]];
    swap_ay[i] = ay[isort[i]];
    swap_zz[i] = zz[isort[i]];
    swap_vz[i] = vz[isort[i]];
    swap_az[i] = az[isort[i]];
    swap_p[i] = p[isort[i]];
    swap_den[i] = den[isort[i]];
    swap_hh[i] = hh[isort[i]];
    swap_mass[i] = mass[isort[i]];
    swap_u[i] = u[isort[i]];
    swap_au[i] = au[isort[i]];
}
for(i = 0; i < n; i++) {
    vx[i] = swap_vx[i];
    ax[i] = swap_ax[i];
    yy[i] = swap_yy[i];
    vy[i] = swap_vy[i];
    ay[i] = swap_ay[i];
    zz[i] = swap_zz[i];
    vz[i] = swap_vz[i];
    az[i] = swap_az[i];
    p[i] = swap_p[i];
    den[i] = swap_den[i];
    hh[i] = swap_hh[i];
    mass[i] = swap_mass[i];
    u[i] = swap_u[i];
    au[i] = swap_au[i];
}

/*********************************************************************
 * Begin main loop.
 ********************************************************************/
// Do until the time (t) exceeds the ending time (tend).
if (debug > 0) printf("beginning main loop\n");
while(t <= tend) {
  if (debug > 0) printf("entering main loop at time = %f\n",t);
    // Store the values of xx, vx, u, and hh at the beginning of the
    // time step to be used later in the particle pushing routine
    // push2. Do the same for y- and z-values.
    for(int i = 0; i < n; i++) {
	xo[i] = xx[i];
	vxo[i] = vx[i];
	yo[i] = yy[i];
	vyo[i] = vy[i];
	zo[i] = zz[i];
	vzo[i] = vz[i];
	uo[i] = u[i];
	ho[i] = hh[i];
	deno[i] = den[i];
    }

    // Pushing the particles 1/2 a time step ahead of the
    // accelerators.
    /* Update values of hh according to Steinmetz and Muller. */
    SM_h(xx, yy, zz, hh, mass, den, divv_x, divv_y, divv_z, t, delt,
	 n);
    if (t != 0) {
	/* Push the particles 1/2 a time step ahead. */
	push1(xx, yy, zz, u, vx, vy, vz, ax, ay, az, au, c1, c2, n);
	// enz_h1(hh, divv, c1, n); /* Advance hh 1/2 a time step
	//			    ahead according to Benz's alg. */
    }

    // Sort on x after pushing.
    r4_sort(xx, isort, partno, n);
    for(i = 0; i < n; i++) {
	swap_vx[i] = vx[isort[i]];
	swap_ax[i] = ax[isort[i]];
	swap_yy[i] = yy[isort[i]];
	swap_vy[i] = vy[isort[i]];
	swap_ay[i] = ay[isort[i]];
	swap_zz[i] = zz[isort[i]];
	swap_vz[i] = vz[isort[i]];
	swap_az[i] = az[isort[i]];
	swap_p[i] = p[isort[i]];
	swap_den[i] = den[isort[i]];
	swap_hh[i] = hh[isort[i]];
	swap_mass[i] = mass[isort[i]];
	swap_u[i] = u[isort[i]];
	swap_au[i] = au[isort[i]];
    }
    for(i = 0; i < n; i++) {
	vx[i] = swap_vx[i];
	ax[i] = swap_ax[i];
	yy[i] = swap_yy[i];
	vy[i] = swap_vy[i];
	ay[i] = swap_ay[i];
	zz[i] = swap_zz[i];
	vz[i] = swap_vz[i];
	az[i] = swap_az[i];
	p[i] = swap_p[i];
	den[i] = swap_den[i];
	hh[i] = swap_hh[i];
	mass[i] = swap_mass[i];
	u[i] = swap_u[i];
	au[i] = swap_au[i];
    }
    
    // Set boundary conditions. The subroutine reflecting creates
    // reflecting (solid wall) boundary conditions.
    if (irefl == 1) {
	reflecting(xx, vx, xo, vxo, yy, vy, yo, vyo, zz, vz, zo, vzo,
		   den, mass, p, u, uo, hh, ho, n);
	/*
	for(int i = 0; i < n; i++) {
	    if (partno[i] <= 50 || partno[i] > n - 50) {
		vx[i] = 0;
		ax[i] = 0;
	    }
	    if (partno[i] <= 50) {
		p[i] = p[50];
		u[i] = u[50];
		den[i] = den[50];
	    } else if (partno[i] > n - 50) {
		p[i] = p[n-49];
		u[i] = u[n-49];
		den[i] = den[n-49];
	    }
	}
	*/
    }

    // The total energy.
    toteo = 0;
    if (t == 0) {
	for(int i = 0; i < n; i++) {
	    vt[i] = 0.5 * sqrt(pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
	}
	if (ientro == 0) {
	    if (irefl == 1) {
		for(int i = 0; i < n; i++) {
		    if (partno[i] > 51 && partno[i] < n - 50) {
			toteo += mass[i] * u[i];
		    }
		}
	    } else {
		for(int i = 0; i < n; i++) {
		    toteo += mass[i] * u[i];
		}
	    }
	} else {
	    if (irefl == 1) {
		for(int i = 0; i < n; i++) {
		    if (partno[i] > 51 && partno[i] < n - 50) {
			toteo += mass[i] * u[i] * pow(den[i],(gamma -
							      1.0)) /
			    (gamma - 1.0);
		    }
		}
	    } else {
		for(int i = 0; i < n; i++) {
		    toteo += mass[i] * u[i] * pow(den[i],(gamma - 1.0)) /
			(gamma - 1.0);
		}
	    }
	}

	if (irefl == 1) {
	    for(int i = 0; i < n; i++) {
		if (partno[i] > 51 && partno[i] < n - 50) {
		    toteo += mass[i] * u[i] * pow(den[i],(gamma - 1.0)) /
			(gamma - 1.0);
		}
	    }
	} else {
	    for(int i = 0; i < n; i++) {
		toteo += mass[i] * vt[i];
	    }
	}
    }

    // Sound speed.
    for(i = 0; i < n; i++) {
	sound[i] = sqrt(gamma * p[i] / den[i]);
    }

    // Thermal diffusion.
    //
    // Here, the subroutine "therm_diffusion" is used, which computes
    // an additional term to be added to the energy equation for each
    // particle as described by Monaghan (1988), unpublished preprint.
    if (itherm == 1) {
	therm_diffusion(xx, vx, yy, vy, zz, vz, hh, mass, den, p, u,
			divv_x, divv_y, divv_z, sound, n, g1, g2,
			gamma, thermdiff);
    }

    // Acclerations and du/dt.
    //
    // The subroutine "accelerations" computes the acceleration of
    // each particle (ax) due to pressure forces and artificial
    // viscosity as well as heating of the particles (au = du/dt).
    accelerations(xx, vx, ax, yy, vy, ay, zz, vz, az, au, hh, p, u,
		  mass, den, divv_x, divv_y, divv_z, sound, maxmu_x,
		  maxmu_y, maxmu_z, n, inorm, alpha, beta, ibulk,
		  alphab, betab, gamma, ientro, ientro2);
    //
    // The subroutine "bouyancy" computes the bouyant force of each
    // particle and adjusts the vertical acceleration (az) so that
    // particles gravitate towards areas of the same density.
    //    bouyancy(xx, yy, zz, az, den, n, b1, grav); //TOBIN

    /*
    if (irefl == 1) {
	for(int i = 0; i < n; i++0) {
	    if (partno[i] >= 50 || partno[i] > n - 50) {
		ax[i] = 0;
	    }
	}
    }
    */

    if (ientro == 1) {
	for(int i = 0; i < n; i++) {
	    au[i] *= u[i] * den[i] * (gamma - 1.0) / p[i];
	}
    }

    // If the switch is set, add the thermal diffusion term.
    if (itherm == 1) {
	for(int i = 0; i < n; i++) {
	    au[i] += thermdiff[i];
	}
    }

    // Push the particles the rest of the way through the time step.
    push2(xx, u, vx, ax, au, c1, c2, delt, xo, vxo, uo, t, n);

    // Update hh, following Benz (1990).
    // benz_h2(hh, delt, divv, ho, n);

    // Sort into x order after pushing particles.
    r4_sort(xx, isort, partno, n);
    for(i = 0; i < n; i++) {
	swap_vx[i] = vx[isort[i]];
	swap_ax[i] = ax[isort[i]];
	swap_p[i] = p[isort[i]];
	swap_den[i] = den[isort[i]];
	swap_hh[i] = hh[isort[i]];
	swap_mass[i] = mass[isort[i]];
	swap_u[i] = u[isort[i]];
	swap_au[i] = au[isort[i]];
    }
    for(i = 0; i < n; i++) {
	vx[i] = swap_vx[i];
	ax[i] = swap_ax[i];
	p[i] = swap_p[i];
	den[i] = swap_den[i];
	hh[i] = swap_hh[i];
	mass[i] = swap_mass[i];
	u[i] = swap_u[i];
	au[i] = swap_au[i];
    }

    // Update the density and pressure with the new positions and new
    // values of the energy density.
    update_den(xx, yy, zz, mass, den, hh, n);
    if (ientro == 1) {
	for(int i = 0; i < n; i++) {
	    p[i] = u[i];
	}
	    
	// If the entropic formulation is used, compute the pressure
	// correspondingly.
	smooth(xx, yy, zz, hh, mass, den, n, p);

	for(i = 0; i < n; i++) {
	    p[i] *= pow(den[i],gamma);
	}
    } else {
	update_p(hh, mass, u, den, p, gamma, n, ientro);
    }
    for(i = 0; i < n; i++) {
	if (p[i] <= 0) p[i] = 0;

	// A new sound speed is computed for each particle.
	sound[i] = sqrt(gamma * p[i] / den[i]);
    }

    // Total energy.
    for(i = 0; i < n; i++) {
	vt[i] = 0.5 * sqrt(pow(vx[i],2) + pow(vy[i],2) + pow(vz[i],2));
    }

    tote = 0;
    if (ientro == 0) {
	if (irefl == 1) {
	    for(int i = 0; i < n; i++) {
		if (partno[i] > 51 && partno[i] < n - 50) {
		    tote += mass[i] * u[i];
		}
	    }
	} else {
	    for(int i = 0; i < n; i++) {
		tote += mass[i] * u[i];
	    }
	}
    } else {
	if (irefl == 1) {
	    for(int i = 0; i < n; i++) {
		if (partno[i] > 51 && partno[i] < n - 50) {
		    tote += mass[i] * u[i] * pow(den[i],(gamma -
							 1.0)) /
			(gamma - 1.0);
		}
	    }
	} else {
	    for(int i = 0; i < n; i++) {
		tote += mass[i] * u[i] * pow(den[i],(gamma - 1.0)) /
		    (gamma - 1.0);
	    }
	}
    }
    
    if (irefl == 1) {
	for(int i = 0; i < n; i++) {
	    if (partno[i] > 51 && partno[i] < n - 50) {
		tote += mass[i] * vt[i];
	    }
	}
    } else {
	for(int i = 0; i < n; i++) {
	    tote += mass[i] * vt[i];
	}
    }
    
    float perc = toteo - tote;
    perc = 100.0 * ((perc < 0) ? -perc : perc) / toteo;
    
    // cout << "tote = " << tote << " perc = " << perc;
    
    // Boundary conditions are reset since the particles have been
    // pushed a second time.
    if (irefl == 1) {
	reflecting(xx, vx, xo, vxo, yy, vy, yo, vyo, zz, vz, zo, vzo, den, mass, p, u, uo, hh, ho, n);
	/*
	for(int i = 0; i < n; i++) {
	    if (partno[i] <= 50 || partno[i] > n - 50) {
		vx[i] = 0;
		ax[i] = 0;
	    }
	    if (partno[i] <= 50) {
		p[i] = p[50];
		u[i] = u[50];
		den[i] = den[50];
	    } else if (partno[i] > n - 50) {
		p[i] = p[n-49];
		u[i] = u[n-49];
		den[i] = den[n-49];
	    }
	    if (partno[i] > n - 50) {
		p[i] = p[n-50];
	    }
	}
	*/
    }

    // Compute the divergence of the velocity field following Monaghan
    // (1988) unpublished preprint.
    div_v(xx, vx, yy, vy, zz, vz, hh, mass, den, n, divv_x, divv_y, divv_z);

    // A new time step is computed according to the courant condition
    // following Hernquist and Katz (1989), also Monaghan (1988).
    new_delt(hh, divv_x, divv_y, divv_z, sound, inorm, alpha, beta,
	     ibulk, alphab, betab, maxmu_x, maxmu_y, maxmu_z, n, c1,
	     c2, cour, delt);

    // OUTPUT

    printf("%f\n",xx[5]);

    // Increment the time.
    t += delt;
}
}
