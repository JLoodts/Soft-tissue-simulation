/*********************************************************************
 * Translated from sph_mods.f in Fortran to C++ by Wesley Miaw
 * <wesley@wesman.net>. Original code by Kevin Olson
 * <olson@jeans.gsfc.nasa.gov> is available from NASA/Goddard at
 * http://sdcd.gsfc.nasa.gov/ESS/exchange/contrib/olson/sph_1d.html.
 *
 * Modified: 05/08/2001
 ********************************************************************/
#include "sph_mods.h"
#include <stdio.h>

static int DEBUG_MODS = 0;

/*********************************************************************
 * The subroutine accelerations computes the acllerations on each
 * particle as well as the rate of change of the specific energy (or
 * entropy), u.
 ********************************************************************/
void accelerations(float *xx, float *vx, float *ax, float *yy, float
		   *vy, float *ay, float *zz, float *vz, float *az,
		   float *au, float *hh, float *p, float *u, float
		   *mass, float *den, float *divv_x, float *divv_y,
		   float *divv_z, float *sound, float *maxmu_x, float
		   *maxmu_y, float *maxmu_z, int n, int inorm, float
		   alpha, float beta, int ibulk, float alphab, float
		   betab, float gamma, int ientro, int ientro2)
{
    if (DEBUG_MODS > 0) printf("acclerations\n");

    // Create some local variables.
    float *xt = new float[n], *vxt = new float[n], *axt = new float[n];
    float *yt = new float[n], *vyt = new float[n], *ayt = new float[n];
    float *zt = new float[n], *vzt = new float[n], *azt = new float[n];
    float *aut = new float[n], *ht = new float[n], *pt = new float[n], *ut = new float[n];
    float *masst = new float[n], *dent = new float[n], *soundt = new float[n];
    float *divv_xt = new float[n], *divv_yt = new float[n], *divv_zt = new float[n];
    float *maxmu_xt = new float[n], *maxmu_yt = new float[n], *maxmu_zt = new float[n];
    float *r = new float[n], *v_x = new float[n], *v_y = new float[n], *v_z = new float[n];
    float *h = new float[n], *gradw = new float[n], *cij = new float[n], *denij = new float[n];
    float *qi_x = new float[n], *qj_x = new float[n], *qi_y = new float[n], *qj_y = new float[n], *qi_z = new float[n], *qj_z = new float[n];
    float *artv_x = new float[n], *artv_y = new float[n], *artv_z = new float[n], *eta = new float[n];
    float *mu_x = new float[n], *mu_y = new float[n], *mu_z = new float[n];
    float *pforce_x = new float[n], *pforce_y = new float[n], *pforce_z = new float[n];

    // Initialize accelerations and du/dt to zero.
    for(int i = 0; i < n; i++) {
	ax[i] = 0;
	ay[i] = 0;
	az[i] = 0;
	au[i] = 0;
    }

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(i = 0; i < n; i++) {
	xt[i] = xx[i];
	vxt[i] = vx[i];
	axt[i] = ax[i];
	yt[i] = yy[i];
	vyt[i] = vy[i];
	ayt[i] = ay[i];
	zt[i] = zz[i];
	vzt[i] = vz[i];
	azt[i] = az[i];
	ht[i] = hh[i];
	aut[i] = au[i];
	pt[i] = p[i];
	ut[i] = u[i];
	masst[i] = mass[i];
	dent[i] = den[i];
	divv_xt[i] = divv_x[i];
	divv_yt[i] = divv_y[i];
	divv_zt[i] = divv_z[i];
	soundt[i] = sound[i];
	maxmu_x[i] = 0;
	maxmu_xt[i] = 0;
	maxmu_y[i] = 0;
	maxmu_yt[i] = 0;
	maxmu_z[i] = 0;
	maxmu_zt[i] = 0;
	r[i] = 0;
    }

    // Compute forces while any two particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Shift the data one location.
	lshift(xt, n);
	lshift(vxt, n);
	lshift(axt, n);
	lshift(yt, n);
	lshift(vyt, n);
	lshift(ayt, n);
	lshift(zt, n);
	lshift(vzt, n);
	lshift(azt, n);
	lshift(ht, n);
	lshift(pt, n);
	lshift(ut, n);
	lshift(aut, n);
	lshift(masst, n);
	lshift(dent, n);
	lshift(divv_xt, n);
	lshift(divv_yt, n);
	lshift(divv_zt, n);
	lshift(soundt, n);
	lshift(maxmu_xt, n);
	lshift(maxmu_yt, n);
	lshift(maxmu_zt, n);

	// Calculate the distance between particles (r), their average
	// smoothing length (h), and their relative velocity (v).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	    h[i] = (ht[i] + hh[i]) / 2.0;
	    v_x[i] = vx[i] - vxt[i];
	    v_y[i] = vy[i] - vyt[i];
	    v_z[i] = vz[i] - vzt[i];
	}

	// Process molecular viscosity.
	if (inorm == 1) {
	    for(int i = 0; i < n; i++) {
		eta[i] = 0.1 * pow(h[i],2);
		mu_x[i] = h[i] * v_x[i] *r[i] / (pow(r[i],2) + eta[i]);
		mu_y[i] = h[i] * v_y[i] *r[i] / (pow(r[i],2) + eta[i]);
		mu_z[i] = h[i] * v_z[i] *r[i] / (pow(r[i],2) + eta[i]);
		if (v_x[i] * r[i] > 0) { mu_x[i] = 0; }
		if (v_y[i] * r[i] > 0) { mu_y[i] = 0; }
		if (v_z[i] * r[i] > 0) { mu_z[i] = 0; }
		if (abs(mu_x[i]) > maxmu_x[i]) { maxmu_x[i] = abs(mu_x[i]); }
		if (abs(mu_y[i]) > maxmu_y[i]) { maxmu_y[i] = abs(mu_y[i]); }
		if (abs(mu_z[i]) > maxmu_z[i]) { maxmu_z[i] = abs(mu_z[i]); }
		cij[i] = (sound[i] + soundt[i]) / 2.0;
		denij[i] = (den[i] + dent[i]) / 2.0;
		artv_x[i] = -alpha * mu_x[i] * cij[i] + beta *
		    pow(mu_x[i],2);
		artv_y[i] = -alpha * mu_y[i] * cij[i] + beta *
		    pow(mu_y[i],2);
		artv_z[i] = -alpha * mu_z[i] * cij[i] + beta *
		    pow(mu_z[i],2);
		artv_x[i] /= denij[i];
		artv_y[i] /= denij[i];
		artv_z[i] /= denij[i];
	    }
	}

	// Process bulk viscosity.
	if (ibulk == 1) {
	    for(int i = 0; i < n; i++) {
		if (divv_x[i] >= 0) {
		    qi_x[i] = 0;
		} else {
		    qi_x[i] = alphab * hh[i] * den[i] * sound[i] *
			abs(divv_x[i]) + betab * pow(hh[i],2) * den[i] *
			pow(divv_x[i],2);
		}
		if (divv_y[i] >= 0) {
		    qi_y[i] = 0;
		} else {
		    qi_y[i] = alphab * hh[i] * den[i] * sound[i] *
			abs(divv_y[i]) + betab * pow(hh[i],2) * den[i] *
			pow(divv_y[i],2);
		}
		if (divv_z[i] >= 0) {
		    qi_z[i] = 0;
		} else {
		    qi_z[i] = alphab * hh[i] * den[i] * sound[i] *
			abs(divv_z[i]) + betab * pow(hh[i],2) * den[i] *
			pow(divv_z[i],2);
		}

		if (divv_xt[i] >= 0) {
		    qj_x[i] = 0;
		} else {
		    qj_x[i] = alphab * ht[i] * dent[i] * soundt[i] *
			abs(divv_xt[i]) + betab * pow(ht[i],2) * dent[i] *
			pow(divv_xt[i],2);
		}
		if (divv_yt[i] >= 0) {
		    qj_y[i] = 0;
		} else {
		    qj_y[i] = alphab * ht[i] * dent[i] * soundt[i] *
			abs(divv_yt[i]) + betab * pow(ht[i],2) * dent[i] *
			pow(divv_yt[i],2);
		}
		if (divv_yt[i] >= 0) {
		    qj_y[i] = 0;
		} else {
		    qj_y[i] = alphab * ht[i] * dent[i] * soundt[i] *
			abs(divv_yt[i]) + betab * pow(ht[i],2) * dent[i] *
			pow(divv_yt[i],2);
		}

		artv_x[i] = qi_x[i] / pow(den[i],2) + qj_x[i] /
		    pow(dent[i],2);
		artv_y[i] = qi_y[i] / pow(den[i],2) + qj_y[i] /
		    pow(dent[i],2);
		artv_z[i] = qi_z[i] / pow(den[i],2) + qj_z[i] /
		    pow(dent[i],2);

		if (r[i] * v_x[i] > 0) { artv_x[i] = 0; }
		if (r[i] * v_y[i] > 0) { artv_y[i] = 0; }
		if (r[i] * v_z[i] > 0) { artv_z[i] = 0; }
	    }
	}

	// Compute the gradient of w.
	gradww(r, hh, ht, n, gradw);

	// Process the forces.
	if (ientro2 == 0) {
	    for(int i = 0; i < n; i++) {
		// The two symmetric forms for the pressure force.
		// pforce[i] = p[i] / den[i]**2 + pt[i] / dent[i]**2;
		pforce_x[i] = 2.0 * sqrt(p[i] * pt[i]) / (den[i] *
							  dent[i]);
		pforce_y[i] = 2.0 * sqrt(p[i] * pt[i]) / (den[i] *
							  dent[i]);
		pforce_z[i] = 2.0 * sqrt(p[i] * pt[i]) / (den[i] *
							  dent[i]);
		
		// Add the art. visc. term and multiply by the gradient of
		// the int. kernal.
		pforce_x[i] = (pforce_x[i] + artv_x[i]) * gradw[i];
		pforce_y[i] = (pforce_y[i] + artv_y[i]) * gradw[i];
		pforce_z[i] = (pforce_z[i] + artv_z[i]) * gradw[i];
		
		// Add to the running totals.
		ax[i] -= masst[i] * pforce_x[i];
		axt[i] += mass[i] * pforce_x[i];
		ay[i] -= masst[i] * pforce_y[i];
		ayt[i] += mass[i] * pforce_y[i];
		az[i] -= masst[i] * pforce_z[i];
		azt[i] += mass[i] * pforce_z[i];
	    }
	} else {
	    for(int i = 0; i < n; i++) {
		// The following is for the case when the entropic
		// formulation of the force is used.
		pforce_x[i] = ut[i] + (gamma - 1.0) * u[i];
		pforce_x[i] *= pow(den[i],(gamma - 2.0));
		pforce_y[i] = pforce_x[i];
		pforce_z[i] = pforce_z[i];
		
		ax[i] -= masst[i] * (pforce_x[i] + artv_x[i]) *
		    gradw[i];
		ay[i] -= masst[i] * (pforce_y[i] + artv_y[i]) *
		    gradw[i];
		az[i] -= masst[i] * (pforce_z[i] + artv_z[i]) *
		    gradw[i];
		pforce_x[i] *= pow(dent[i],(gamma - 2.0));
		pforce_y[i] *= pow(dent[i],(gamma - 2.0));
		pforce_z[i] *= pow(dent[i],(gamma - 2.0));
		
		axt[i] += mass[i] * (pforce_x[i] + artv_x[i]) *
		    gradw[i];
		ayt[i] += mass[i] * (pforce_y[i] + artv_y[i]) *
		    gradw[i];
		azt[i] += mass[i] * (pforce_z[i] + artv_z[i]) *
		    gradw[i];
	    }
	}

	// Process du/dt.
	if (ientro == 0) {
	    for(int i = 0; i < n; i++) {
		au[i] += masst[i] * (pforce_x[i] * v_x[i] +
				     pforce_y[i] * v_y[i] +
				     pforce_z[i] * v_z[i]);
		aut[i] += mass[i] * (pforce_x[i] * v_x[i] +
				     pforce_y[i] * v_y[i] +
				     pforce_z[i] * v_z[i]);
	    }
	} else {
	    for(int i = 0; i < n; i++) {
		// The rate of change of the entropic function.
		au[i] += masst[i] * (artv_x[i] * v_x[i] + artv_y[i] *
				     v_y[i] + artv_z[i] * v_z[i]) *
		    gradw[i];
		aut[i] += mass[i] * (artv_x[i] * v_x[i] + artv_y[i] *
				     v_y[i] + artv_z[i] * v_z[i]) *
		    gradw[i];
	    }
	}

	// Increment the counter.
	count++;
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(axt, n);
	rshift(ayt, n);
	rshift(azt, n);
	rshift(aut, n);
	rshift(maxmu_xt, n);
	rshift(maxmu_yt, n);
	rshift(maxmu_zt, n);
    }
    for(i = 0; i < n; i++) {
	ax[i] += axt[i];
	ay[i] += ayt[i];
	az[i] += azt[i];
	au[i] = 0.5 * (au[i] + aut[i]);
	if (maxmu_xt[i] > maxmu_x[i]) { maxmu_x[i] = maxmu_xt[i]; }
	if (maxmu_yt[i] > maxmu_y[i]) { maxmu_y[i] = maxmu_yt[i]; }
	if (maxmu_zt[i] > maxmu_z[i]) { maxmu_z[i] = maxmu_zt[i]; }
    }
}

/*********************************************************************
 * The subroutine bouyancy adjusts the vertical accleration (az) of
 * each particle so that particles gravitate towards other particles
 * of the same density. Less dense particles should gravitate upwards
 * and denser particles should gravitate downwards.
 *
 * Assumes the x positions are sorted in increasing order.
 ********************************************************************/
void bouyancy(float *xx, float *yy, float *zz, float *az, float *den,
	      int n, float b1, float grav)
{
    if (DEBUG_MODS > 0) printf("bouyancy\n");

    // For each particle 1 through n, compare that particle to all
    // particles with the same x and y coordinate values who is also
    // either one position above or one position below.
    for(int i = 0; i < n; i++) {
	for(int j = i+1; j < n; j++) {
	    // Don't search any farther ahead if the x-values no
	    // match. The x-values should be sorted.
	    if (xx[i] != xx[j]) break;

	    // Otherwise, if the y-values are identical and the
	    // z-values differ by only 1, adjust the particle
	    // acceleration appropriately.
	    if (yy[i] == yy[j] && abs(zz[i] - zz[j]) == 1) {
		// If the density above is higher, rise.
		// Otherwise, if the density below is less, sink.
		if (zz[j] > zz[i] && den[j] > den[i]) {
			az[i] += b1 * (den[j] - den[i]);
		} else if (zz[j] < zz[i] && den[j] < den[i]) {
			az[i] += (b1 + grav) * (den[j] - den[i]);
		}
	    }
	}
    }
}

/*********************************************************************
 * The subroutine therm_diffusion computes an additional thermal
 * diffusion term which can be added to the energy equation.
 ********************************************************************/
void therm_diffusion(float *xx, float *vx, float *yy, float *vy, float
		     *zz, float *vz, float *hh, float *mass, float
		     *den, float *p, float *u, float *divv_x, float
		     *divv_y, float *divv_z, float *sound, int n,
		     float g1, float g2, float gamma, float
		     *thermdiff)
{
    if (DEBUG_MODS > 0) printf("therm_diffusion\n");

    // Create some local variables.
    float *xt = new float[n], *vxt = new float[n], *yt = new float[n], *vyt = new float[n], *zt = new float[n], *vzt = new float[n];
    float *ht = new float[n], *pt = new float[n], *ut = new float[n], *masst = new float[n], *dent = new float[n];
    float *divv_xt = new float[n], *divv_yt = new float[n], *divv_zt = new float[n];
    float *soundt = new float[n], *thermdifft = new float[n];
    float *r = new float[n], *uu = new float[n], *gradw = new float[n], *qqi = new float[n], *qqj = new float[n];
    float *denij = new float[n], *h = new float[n], *eta = new float[n];

    // Initialize the thermal diffusion to zero.
    for(int i = 0; i < n; i++) {
	thermdiff[i] = 0;
	thermdifft[i] = 0;
    }

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(i = 0; i < n; i++) {
	xt[i] = xx[i];
	vxt[i] = vx[i];
	yt[i] = yy[i];
	vyt[i] = vy[i];
	zt[i] = zz[i];
	vzt[i] = vz[i];
	ht[i] = hh[i];
	pt[i] = p[i];
	ut[i] = u[i];
	masst[i] = mass[i];
	dent[i] = den[i];
	divv_xt[i] = divv_x[i];
	divv_yt[i] = divv_y[i];
	divv_zt[i] = divv_z[i];
	sound[i] = sqrt(gamma * p[i] / den[i]);
	soundt[i] = sound[i];
	r[i] = 0;
    }

    // Compute thermal diffusion while any two particles are
    // overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Shift the data one location.
	lshift(xt, n);
	lshift(vxt, n);
	lshift(yt, n);
	lshift(vyt, n);
	lshift(zt, n);
	lshift(vzt, n);
	lshift(ht, n);
	lshift(pt, n);
	lshift(ut, n);
	lshift(masst, n);
	lshift(dent, n);
	lshift(divv_xt, n);
	lshift(divv_yt, n);
	lshift(divv_zt, n);
	lshift(soundt, n);
	lshift(thermdifft, n);

	// Calculate the distance between particles (r), their
	// relative energy (uu), (qqi), (qqj), (denij).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	    uu[i] = u[i] - ut[i];

	    qqi[i] = g1 * hh[i] * den[i] * sound[i] + g2 * den[i] *
		pow(hh[i],2) * ((abs(divv_x[i]) - divv_x[i]) +
				(abs(divv_y[i]) - divv_y[i]) +
				(abs(divv_z[i]) - divv_z[i]));
	    qqi[i] /= den[i];

	    qqj[i] = g1 * ht[i] * dent[i] * soundt[i] + g2 * dent[i] *
		pow(ht[i],2) * ((abs(divv_xt[i]) - divv_xt[i]) +
				(abs(divv_yt[i]) - divv_yt[i]) +
				(abs(divv_zt[i]) - divv_zt[i]));
	    qqj[i] /= dent[i];
	    qqj[i] = (qqi[i] + qqj[i]) / 2.0;

	    denij[i] = (den[i] + dent[i]) / 2.0;
	}

	// Compute the gradient of w.
	gradww(r, hh, ht, n, gradw);

	for(i = 0; i < n; i++) {
	    h[i] = (hh[i] + ht[i]) / 2.0;
	    eta[i] = 0.1 * pow(h[i],2);
	    
	    if (r[i] != 0) {
		thermdiff[i] += masst[i] * qqj[i] * uu[i] * r[i] *
		    gradw[i] / (denij[i] * (pow(r[i],2) + eta[i]));
		thermdifft[i] -= mass[i] * qqj[i] * uu[i] * r[i] *
		    gradw[i] / (denij[i] * (pow(r[i],2) + eta[i]));
	    }
	}

	// Increment the counter.
	count++;
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(thermdifft, n);
    }
    for(i = 0; i < n; i++) {
	thermdiff[i] = 2.0 * (thermdiff[i] + thermdifft[i]);
    }
}

/*********************************************************************
 * The subroutine push1 advances the particles' positions, velocities,
 * and specific energies forward in time by half a time step.
 ********************************************************************/
void push1(float *xx, float *yy, float *zz, float *u, float *vx, float
	   *vy, float *vz, float *ax, float *ay, float *az, float *au,
	   float c1, float c2, int n)
{
    if (DEBUG_MODS > 0) printf("push1\n");

    for(int i = 0; i < n; i++) {
	vx[i] += c1 * ax[i];
	vy[i] += c1 * ay[i];
	vz[i] += c1 * az[i];
	xx[i] += c2 * ax[i] + c1 * vx[i];
	yy[i] += c2 * ay[i] + c1 * vy[i];
	zz[i] += c2 * az[i] + c1 * vz[i];
	u[i] += c1 * au[i];
    }
}

/*********************************************************************
 * The subroutine benz_h1 advances the particles' smoothing length
 * forward in time by half a time step according to the algorithm of
 * Benz.
 ********************************************************************/
void benz_h1(float *hh, float *divv, float c1, int n)
{
    if (DEBUG_MODS > 0) printf("benz_h1\n");

    for(int i = 0; i < n; i++) {
	hh[i] += c1 * hh[i] * divv[i];
    }
}

/*********************************************************************
 * The subroutine push2 is called after the accelerations have been
 * computed and advances the particles through a complete time step.
 ********************************************************************/
void push2(float *xx, float *u, float *vx, float *ax, float *au, float
	   c1, float c2, float delt, float *xo, float *vxo, float *uo,
	   float t, int n)
{
    if (DEBUG_MODS > 0) printf("push2\n");

    if (t == 0) {
	for(int i = 0; i < n; i++) {
	    xx[i] += c1 * vx[i];
	    vx[i] += c1 * ax[i];
	    u[i] += c1 * au[i];
	}
    } else {
	for(int i = 0; i < n; i++) {
	    vx[i] = vxo[i] + delt * ax[i];
	    xx[i] = xo[i] + delt * (vxo[i] + c1 * ax[i]);
	    u[i] = uo[i] + delt * au[i];
	}
    }
}

/*********************************************************************
 * The subroutine benz_h2 advances the smoothing lengths through a
 * complete time step and is called in consort with push2.
 ********************************************************************/
void benz_h2(float *hh, float c1, float delt, float *divv, float *ho,
	     float t, int n)
{
    if (DEBUG_MODS > 0) printf("benz_h2\n");

    if (t == 0) {
	for(int i = 0; i < n; i++) {
	    hh[i] += c1 * divv[i] * hh[i];
	}
    } else {
	for(int i = 0; i < n; i++) {
	    hh[i] = ho[i] + delt * divv[i] * hh[i];
	}
    }
}

/*********************************************************************
 * The subroutine update_den computes the density according to the
 * positions and the smoothing lengths of the particles.
 ********************************************************************/
void update_den(float *xx, float *yy, float *zz, float *mass, float
		*den, float *hh, int n)
{
    if (DEBUG_MODS > 0) printf("update_den\n");

    // Create some local variables.
    float *xt = new float[n], *yt = new float[n], *zt = new float[n], *dent = new float[n], *masst = new float[n], *ht = new float[n];
    float *w = new float[n], *r = new float[n];

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(int i = 0; i < n; i++) {
	xt[i] = xx[i];
	yt[i] = yy[i];
	zt[i] = zz[i];
	ht[i] = hh[i];
	masst[i] = mass[i];
	den[i] = 0;
	dent[i] = 0;
	r[i] = 0;
    }

    // Compute densities while any two particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Calculate the distance between particles (r), (w).
	if (DEBUG_MODS > 2) printf("calculating r[i]\n");
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	}
	ww(r, hh, ht, n, w);

	// Keep a running total of the density (den) and the shifted
	// density (dent).
	for(i = 0; i < n; i++) {
	    den[i] += masst[i] * w[i];
	    if (count > 0) { dent[i] += mass[i] * w[i]; }
	}

	// Increment the counter.
	count++;

	// Shift the data.
	lshift(xt, n);
	lshift(yt, n);
	lshift(zt, n);
	lshift(ht, n);
	lshift(masst, n);
	lshift(dent, n);
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(dent, n);
    }
    for(i = 0; i < n; i++) {
	den[i] += dent[i];
    }
}

/*********************************************************************
 * The subroutine update_p computes the pressure for each particle
 * assuming an ideal gas.
 ********************************************************************/
void update_p(float *hh, float *mass, float *u, float *den, float *p,
	      float gamma, int n, int ientro)
{
    if (DEBUG_MODS > 0) printf("update_p\n");

    if (ientro == 0) {
	for(int i = 0; i < n; i++) {
	    p[i] = (gamma - 1.0) * u[i] * den[i];
	}
    } else {
	for(int i = 0; i < n; i++) {
	    p[i] = u[i] * pow(den[i],gamma);
	}
    }
}

/*********************************************************************
 * The subroutine new_delt computes a new time step.
 ********************************************************************/
void new_delt(float *hh, float *divv_x, float *divv_y, float *divv_z,
	      float *sound, int inorm, float alpha, float beta, int
	      ibulk, float alphab, float betab, float *maxmu_x, float
	      *maxmu_y, float *maxmu_z, int n, float c1, float c2,
	      float cour, float delt)
{
    if (DEBUG_MODS > 0) printf("new_delt\n");

    // Create some local variables.
    float *delti = new float[n], *ht = new float[n];

    // In the case the molecular viscosity is used.
    if (inorm == 1 ||
	(inorm == 0 && ibulk == 0))
    {
	for(int i = 0; i < n; i++) {
	    delti[i] = hh[i] * min(abs(divv_x[i]), min(abs(divv_y[i]),
						       abs(divv_z[i]))) +
		sound[i] + 1.2 * (alpha * sound[i] + beta *
				  min(maxmu_x[i], min(maxmu_y[i],
						      maxmu_z[i])));
	    delti[i] = hh[i] / delti[i];
	}
    }

    // In the case the bulk viscosity is used.
    if (ibulk == 1) {
	for(int i = 0; i < n; i++) {
	    ht[i] = betab * hh[i] * min(abs(divv_x[i]),
					min(abs(divv_y[i]), abs(divv_z[i])));
	    if (min(divv_x[i], min(divv_y[i], divv_z[i])) >= 0) {
		ht[i] = 0;
	    }
	    delti[i] = hh[i] * min(abs(divv_x[i]),
				   min(abs(divv_y[i]),
				       abs(divv_z[i]))) + sound[i] +
		1.2 * (alphab * sound[i] + ht[i]);
	    delti[i] = hh[i] / delti[i];
	}
    }

    delt = minval(delti, n);
    delt *= cour;
    c1 = delt / 2.0;
    c2 = pow(delt,2) / 4.0;
}

/*********************************************************************
 * The subroutine reflecting sets the boundary conditions. This is
 * done by using 50 guard particles on either end of the domain which
 * are the mirror images of particles in the domain.
 ********************************************************************/
void reflecting(float *xx, float *vx, float *xo, float *vxo, float
		*yy, float *vy, float *yo, float *vyo, float *zz,
		float *vz, float *zo, float *vzo, float *den, float
		*mass, float *p, float *u, float *uo, float *hh, float
		*ho, int n)
{
    if (DEBUG_MODS > 0) printf("reflecting\n");

    for(int i = 0; i < 50; i++) {
	xx[i] = -xx[100-i] + xx[50];
	vx[i] = -vx[100-i];
	xo[i] = -xo[100-i] + xo[50];
	vxo[i] = -vxo[100-i];
	yy[i] = -yy[100-i] + yy[50];
	vy[i] = -vy[100-i];
	yo[i] = -yo[100-i] + yo[50];
	vyo[i] = -vyo[100-i];
	zz[i] = -zz[100-i] + zz[50];
	vz[i] = -vz[100-i];
	zo[i] = -zo[100-i] + zo[50];
	vzo[i] = -vzo[100-i];
	den[i] = den[100-i];
	mass[i] = mass[100-i];
	p[i] = p[100-i];
	u[i] = u[100-i];
	uo[i] = uo[100-i];
	hh[i] = hh[100-i];
	ho[i] = ho[100-i];
    }
    int j = n - 53;
    for(i = n - 51; i < n && j > n - 103; i++, j--) {
	xx[i] = xx[n-52] + xx[n-52] - xx[j];
	vx[i] = -vx[j];
	xo[i] = xo[n-52] + xo[n-52] - xo[j];
	vxo[i] = -vxo[j];
	yy[i] = yy[n-52] + yy[n-52] - yy[j];
	vy[i] = -vy[j];
	yo[i] = yo[n-52] + yo[n-52] - yo[j];
	vyo[i] = -vyo[j];
	zz[i] = zz[n-52] + zz[n-52] - zz[j];
	vz[i] = -vz[j];
	zo[i] = zo[n-52] + zo[n-52] - zo[j];
	vzo[i] = -vzo[j];
	den[i] = den[j];
	mass[i] = mass[j];
	p[i] = p[j];
	u[i] = u[j];
	hh[i] = hh[j];
    }
}

/*********************************************************************
 * The subroutine ww computes the 1-d value of the interpolating
 * kernal given a distance (r) between the particles and their
 * smoothing lengths.
 ********************************************************************/
void ww(float *r, float *hh, float *ht, int n, float *w)
{
    if (DEBUG_MODS > 0) printf("ww\n");

    // Create some local variables.
    float q1 = 1.0;
    float q2 = 2.0;
    float *q = new float[n], *w2 = new float[n];

    for(int i = 0; i < n; i++) {
	q[i] = ((r[i] < 0) ? -r[i] : r[i]) / hh[i];
	
	if (q[i] <= q1) {
	    w[i] = 2.0 / 3.0 - pow(q[i],2) + 0.5 * pow(q[i],3);
	    w[i] /= hh[i];
	} else if (q[i] > q1 && q[i] <= q2) {
	    w[i] = pow(2.0 - q[i],3);
	    w[i] /= (6.0 * hh[i]);
	} else if (q[i] > q2) {
	    w[i] = 0;
	}

	q[i] = ((r[i] < 0) ? -r[i] : r[i]) / ht[i];

	if (q[i] <= q1) {
	    w2[i] = 2.0 / 3.0 - pow(q[i],2) + 0.5 * pow(q[i],3);
	    w2[i] /= ht[i];
	} else if (q[i] > q1 && q[i] <= q2) {
	    w2[i] = pow(2.0 - q[i],3);
	    w2[i] /= (6.0 * ht[i]);
	} else if (q[i] > q2) {
	    w2[i] = 0;
	}
	
	w[i] = (w[i] + w2[i]) / 2.0;
    }
}

/*********************************************************************
 * The subroutine gradww computes the gradient of the interpolating
 * kernal.
 ********************************************************************/
void gradww(float *r, float *hh, float *ht, int n, float *gradw)
{
    if (DEBUG_MODS > 0) printf("gradww\n");

    // Create some local variables.
    float  *gradw2 = new float[n], *q = new float[n];
    float q1 = 1.0;
    float q2 = 2.0;

    for(int i = 0; i < n; i++) {
	q[i] = ((r[i] < 0) ? -r[i] : r[i]) / hh[i];

	if (q[i] <= q1) {
	    gradw[i] = 1.5 * pow(r[i],2) / pow(hh[i],3);
	    gradw[i] -= 2.0 * ((r[i] < 0) ? -r[i] : r[i]) / pow(hh[i],2);
	    gradw[i] /= hh[i];
	} else if (q[i] > q1 && q[i] <= q2) {
	    gradw[i] = pow(2.0 - q[i],2);
	    gradw[i] = -3.0 * gradw[i] / hh[i];
	    gradw[i] = gradw[i] / (6.0 * hh[i]);
	}
	
	if (r[i] != 0) {
	    gradw[i] *= r[i] / ((r[i] < 0) ? -r[i] : r[i]);
	} else {
	    gradw[i] = 0;
	}

	if (q[i] > q2) { gradw[i] = 0; }

	q[i] = ((r[i] < 0) ? -r[i] : r[i]) / ht[i];

	if (q[i] <= q1) {
	    gradw2[i] = 1.5 * pow(r[i],2) / pow(ht[i],3);
	    gradw2[i] -= 2.0 * ((r[i] < 0) ? -r[i] : r[i]) / pow(ht[i],2);
	    gradw2[i] /= ht[i];
	} else if (q[i] > q1 && q[i] <= q2) {
	    gradw2[i] = pow(2.0 - q[i],2);
	    gradw2[i] = -3.0 * gradw[i] / ht[i];
	    gradw2[i] = gradw[i] / (6.0 * ht[i]);
	}

	if (r[i] != 0) {
	    gradw2[i] *= r[i] / ((r[i] < 0) ? -r[i] : r[i]);
	} else {
	    gradw2[i] = 0;
	}

	if (q[i] > q2) { gradw2[i] = 0; }

	gradw[i] = (gradw[i] + gradw2[i]) / 2.0;
    }
}

/*********************************************************************
 * The subroutine smooth smooths a variable (here represented by u).
 ********************************************************************/
void smooth(float *xx, float *yy, float *zz, float *hh, float *mass,
	    float *den, int n, float *u)
{
    if (DEBUG_MODS > 0) printf("smooth\n");

    // Create some local variables.
    float *xt = new float[n], *yt = new float[n], *zt = new float[n], *ht = new float[n], *masst = new float[n];
    float *old_u = new float[n], *old_ut = new float[n], *ut = new float[n];
    float *r = new float[n], *w = new float[n];

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(int i = 0; i < n; i++) {
	xt[i] = xx[i];
	yt[i] = yy[i];
	zt[i] = zz[i];
	ht[i] = hh[i];
	masst[i] = mass[i];
	old_u[i] = u[i];
	old_ut[i] = u[i];
	u[i] = 0;
	ut[i] = 0;
	r[i] = 0;
    }

    // Smooth the variable while any two particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Calculate the distance between particles (r).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	}
	ww(r, hh, ht, n, w);

	// Smooth the variable.
	for(i = 0; i < n; i++) {
	    u[i] += masst[i] * old_ut[i] * w[i];
	}
	if (count > 0) {
	    for(int i = 0; i < n; i++) {
		ut[i] += mass[i] * old_u[i] * w[i];
	    }
	}
	
	// Increment the counter.
	count++;

	// Shift the data one location.
	lshift(xt, n);
	lshift(yt, n);
	lshift(zt, n);
	lshift(ht, n);
	lshift(masst, n);
	lshift(old_ut, n);
	lshift(ut, n);
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(ut, n);
    }
    for(i = 0; i < n; i++) {
	u[i] += ut[i];
	u[i] /= den[i];
    }
}

/*********************************************************************
 * The subroutine div_v computes the divergence of the velocity field.
 ********************************************************************/
void div_v(float *xx, float *vx, float *yy, float *vy, float *zz,
	   float *vz, float *hh, float *mass, float *den, int n,
	   float *divv_x, float *divv_y, float *divv_z)
{
    if (DEBUG_MODS > 0) printf("div_v\n");

    // Create local variables.
    float *xt = new float[n], *vxt = new float[n], *yt = new float[n], *vyt = new float[n], *zt = new float[n], *vzt = new float[n];
    float *ht = new float[n], *divv_xt = new float[n], *divv_yt = new float[n], *divv_zt = new float[n], *masst = new float[n];
    float *r = new float[n], *gradw = new float[n], *v_x = new float[n], *v_y = new float[n], *v_z = new float[n];

    // Initialize the divergence to zero.
    for(int i = 0; i < n; i++) {
	divv_x[i] = 0;
	divv_y[i] = 0;
	divv_z[i] = 0;
    }

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(i = 0; i < n; i++) {
	xt[i] = xx[i];
	vxt[i] = vx[i];
	yt[i] = yy[i];
	vyt[i] = vy[i];
	zt[i] = zz[i];
	vzt[i] = vz[i];
	ht[i] = hh[i];
	divv_xt[i] = 0;
	divv_yt[i] = 0;
	divv_zt[i] = 0;
	masst[i] = mass[i];
	r[i] = 0;
    }

    // Compute divergence while any two particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Shift the data one location.
	lshift(xt, n);
	lshift(vxt, n);
	lshift(yt, n);
	lshift(vyt, n);
	lshift(zt, n);
	lshift(vzt, n);
	lshift(ht, n);
	lshift(divv_xt, n);
	lshift(divv_yt, n);
	lshift(divv_zt, n);
	lshift(masst, n);

	// Calculate the distance between particles (r), and their
	// relative velocity (v).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	    v_x[i] = vx[i] - vxt[i];
	    v_y[i] = vy[i] - vyt[i];
	    v_z[i] = vz[i] - vzt[i];
	}

	// Compute the gradient of w.
	gradww(r, hh, ht, n, gradw);

	// Compute the divergence. v is in the opposite sense as well
	// as gradw, hence the +.
	for(i = 0; i < n; i++) {
	    divv_x[i] += masst[i] * v_x[i] * gradw[i];
	    divv_xt[i] += mass[i] * v_x[i] * gradw[i];
	    divv_y[i] += masst[i] * v_y[i] * gradw[i];
	    divv_yt[i] += mass[i] * v_y[i] * gradw[i];
	    divv_z[i] += masst[i] * v_z[i] * gradw[i];
	    divv_zt[i] += mass[i] * v_z[i] * gradw[i];
	}

	// Increment the counter.
	count++;
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(divv_xt, n);
	rshift(divv_yt, n);
	rshift(divv_zt, n);
    }
    for(i = 0; i < n; i++) {
	divv_x[i] = (divv_x[i] + divv_xt[i]) / den[i];
	divv_y[i] = (divv_y[i] + divv_yt[i]) / den[i];
	divv_z[i] = (divv_z[i] + divv_zt[i]) / den[i];
    }
}

/*********************************************************************
 * The subroutine SM_h computes new smoothing lengths according to the
 * alg. of Steinmetx and Muller.
 ********************************************************************/
void SM_h(float *xx, float *yy, float *zz, float *hh, float *mass,
	  float *den, float *divv_x, float *divv_y, float *divv_z,
	  float t, float delt, int n)
{
    if (DEBUG_MODS > 0) printf("SM_h\n");

    // Create some local variables.
    float *densmooth = new float[n], *dsmoothdt = new float[n], *avh = new float[n];
    float aveden;

    // Compute a smoothed value of the density.
    smooth_den(xx, yy, zz, hh, mass, den, n, densmooth);

    // Compute the rate of change of the smoothed density.
    d_smooth_den(xx, yy, zz, hh, mass, den, densmooth, divv_x, divv_y,
		 divv_z, n, dsmoothdt);

    // Estimate the smoothed density.
    if (t == 0) {
	for(int i = 0; i < n; i++) {
	    densmooth[i] = den[i];
	}
    }
    for(int i = 0; i < n; i++) {
	densmooth[i] += delt * dsmoothdt[i];
    }

    // Calculate the average density of the system.
    for(i = 0; i < n; i++) {
	densmooth[i] = 1.0 / densmooth[i];
    }
    aveden = 0;
    for(i = 0; i < n; i++) {
	aveden += densmooth[i];
    }
    aveden = 1.0 / aveden;

    // Find new smoothing length.
    float hfac = 1.5;
    for(i = 0; i < n; i++) {
	avh[i] = hfac * mass[i] / aveden;
    }

    // For 3-d need to use (mass/den)**1/3 and (mass/den)**1/2 for
    // 2-d.
    for(i = 0; i < n; i++) {
	densmooth[i] = 1.0 / densmooth[i];
	float kkk = 3.0;
	hh[i] = avh[i] * pow(aveden / densmooth[i],kkk);
    }
}

/*********************************************************************
 * The subroutine smooth_den computes a smoothed value for the density
 * for the Steinmetz and Muller alg.
 ********************************************************************/
void smooth_den(float *xx, float *yy, float *zz, float *hh, float
		*mass, float *den, int n, float *densmooth)
{
    if (DEBUG_MODS > 0) printf("smooth_den\n");

    // Create some local variables.
    float *xt = new float[n], *yt = new float[n], *zt = new float[n], *ht = new float[n], *masst = new float[n], *dent = new float[n], *denst = new float[n];
    float *r = new float[n], *w = new float[n];
    float hfac = 1.5;

    // A counter.
    int count = 0;

    // Store copies of paramters for circular shift.
    for(int i = 0; i < n; i++) {
	xt[i] = xx[i];
	yt[i] = yy[i];
	zt[i] = zz[i];
	ht[i] = hh[i];
	masst[i] = mass[i];
	dent[i] = den[i];
	densmooth[i] = 0;
	denst[i] = 0;
	r[i] = 0;
    }

    // Compute the smoothed value for the density while any two
    // particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Calculate the relative distance between particles (r).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt(pow(xx[i] - xt[i],2) +
			pow(yy[i] - yt[i],2) +
			pow(zz[i] - zt[i],2));
	}
	ww(r, hh, ht, n, w);

	// Compute the smoothed value.
	for(i = 0; i < n; i++) {
	    densmooth[i] += masst[i] * w[i] * dent[i];
	    if (count > 0) { denst[i] += mass[i] * w[i] * den[i]; }
	}

	// Increment the counter.
	count++;

	// Shift the data one location.
	lshift(xt, n);
	lshift(yt, n);
	lshift(zt, n);
	lshift(ht, n);
	lshift(masst, n);
	lshift(dent, n);
	lshift(denst, n);
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < n; i++) {
	rshift(denst, n);
    }
    for(i = 0; i < n; i++) {
	densmooth[i] = (densmooth[i] + denst[i]) / den[i];
	hh[i] = hfac * mass[i] / densmooth[i];
    }
}

/*********************************************************************
 * The subroutine d_smooth_den computes the rate of change of the
 * smoothed density for the Steinmetz + Muller alg.
 ********************************************************************/
void d_smooth_den(float *xx, float *yy, float *zz, float *hh, float
		  *mass, float *den, float *densmooth, float *divv_x,
		  float *divv_y, float *divv_z, int n, float
		  *dsmoothdt)
{
    if (DEBUG_MODS > 0) printf("d_smooth_den\n");

    // Create some local variables.
    float *xt = new float[n], *yt = new float[n], *zt = new float[n];
    float *ht = new float[n], *masst = new float[n], *dent = new float[n], *dsmt = new float[n];
    float *divv_xt = new float[n], *divv_yt = new float[n], *divv_zt = new float[n];
    float *w = new float[n], *r = new float[n];

    // A counter.
    int count = 0;

    // Store copies of parameters for circular shift.
    for(int i = 0; i < n; i++) {
	xt[i] = xx[i];
	yt[i] = yy[i];
	zt[i] = zz[i];
	ht[i] = hh[i];
	masst[i] = mass[i];
	dent[i] = densmooth[i];
	divv_xt[i] = divv_x[i];
	divv_yt[i] = divv_y[i];
	divv_zt[i] = divv_z[i];
	dsmoothdt[i] = 0;
	dsmt[i] = 0;
	r[i] = 0;
    }

    // Compute the rate of change of the smoothed density while any
    // two particles are overlapping.
    while(any_overlap(r, hh, ht, n)) {
	// Calculate the distance between particles (r).
	for(int i = 0; i < n; i++) {
	    r[i] = sqrt((xx[i] - xt[i]) +
			(yy[i] - yt[i]) +
			(zz[i] - zt[i]));
	}
	ww(r, hh, ht, n, w);

	// Compute the rate of change.
	for(i = 0; i < n; i++) {
	    dsmoothdt[i] += masst[i] * dent[i] * (divv_xt[i] +
						  divv_yt[i] +
						  divv_zt[i]) * w[i];
	    if (count > 0) {
		dsmt[i] += mass[i] * densmooth[i] * (divv_x[i] +
						     divv_y[i] +
						     divv_z[i]) *
		    w[i];
	    }
	}

	// Increment the counter.
	count++;

	// Shift the data one location.
	lshift(xt, n);
	lshift(yt, n);
	lshift(zt, n);
	lshift(ht, n);
	lshift(masst, n);
	lshift(dent, n);
	lshift(divv_xt, n);
	lshift(divv_yt, n);
	lshift(divv_zt, n);
	lshift(dsmt, n);
    }

    // Shift the temporary totals back to their starting points and
    // add to get the complete sum.
    for(i = 0; i < count; i++) {
	rshift(dsmt, n);
    }
    for(i = 0; i < count; i++) {
	dsmoothdt[i] += dsmt[i];
	dsmoothdt[i] *= 2.0 / den[i];
	dsmoothdt[i] = densmooth[i] * (divv_x[i] + divv_y[i] +
				       divv_z[i]) - dsmoothdt[i];
    }
}

/*********************************************************************
 * The subroutine any_overlap returns true if any two particles are
 * overlapping.
 ********************************************************************/
bool any_overlap(float *r, float *hh, float *ht, int n) {
    if (DEBUG_MODS > 1) printf("any_overlap(%i)\n", n);

    for(int i = 0; i < n; i++) {
	if (DEBUG_MODS > 2) printf("%i.", i);
	float abs_r = abs(r[i]);
	if (abs_r <= 2*hh[i] || abs_r <= 2*ht[i]) {
	    return true;
	}
    }
    if (DEBUG_MODS > 2) printf("\n");
    return false;
}

/*********************************************************************
 * The subroutine lshift performs a left circular shift of the array
 * by one position.
 ********************************************************************/
void lshift(float *a, int n) {
    if (DEBUG_MODS > 1) printf("lshift\n");

    float old_first = a[0];
    for(int i = 0; i < n-1; i++) {
	a[i] = a[i+1];
    }
    a[n-1] = old_first;
}

/*********************************************************************
 * The subroutine rshift performs a right circular shift of the array
 * by one position.
 ********************************************************************/
void rshift(float *a, int n) {
    if (DEBUG_MODS > 1) printf("rshift\n");

    float old_last = a[n-1];
    for(int i = n-1; i > 0; i--) {
	a[i] = a[i-1];
    }
    a[0] = old_last;
}

/*********************************************************************
 * The subroutine minval returns the minimum value of the array.
 ********************************************************************/
float minval(float *a, int n) {
    if (DEBUG_MODS > 1) printf("minval\n");

    float min = a[0];
    for(int i = 1; i < n; i++) {
	if (a[i] < min) { min = a[i]; }
    }
    return min;
}

/*********************************************************************
 * The subroutine abs returns the absolute value of the float.
 ********************************************************************/
float abs(float x) {
    return (x < 0.0) ? -x : x;
}

/*********************************************************************
 * The subroutine min returns the minimum of two floats.
 ********************************************************************/
float min(float x, float y) {
    return (x < y) ? x : y;
}
