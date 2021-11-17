void noh(double xx[],double vx[],double mass[],double den[],
		 double p[],double u[],double hh[],double gamma,
		 double xmax,double xmin,int partno[],double hfac,int n,
		 bool ientro)
{
	/* This subroutine sets up the initial conditions for the `shocks of Noh',
	 * i.e. a discontinuity in velocity.
	 */

	for(int i=0; i<n; ++i) { 
		p[i]	= 0.01;
		den[i]	= 1.0;
		vx[i]	= 0.0;
		u[i]	= p[i]/den[i];
		u[i]	= u[i]/(gamma-1.0);

		if (ientro==1) u[i] = p[i]/pow(den[i],gamma);

		partno[i] = i;
	}
	double delx	= (xmax-xmin)/n;

	for(i=0; i<n; ++i) { 
		hh[i]	= hfac*delx;
		mass[i]	= den[i]*delx;
	
		if (i==1) {
			xx[i] = 0.0;
		}
		else {
			xx[i] = xx[i-1]+delx;
		}
		if (i<=(n/2)) {
			vx[i] = 1.0;
		}
		else {
			vx[i] = -1.0;
		}
	}
}