void sod_const (double xx[],double vx[],double mass[],
				double den[],double p[],double u[],double hh[],
				double gamma,double xmax,double xmin,int partno[],
				double hfac,int n,bool ientro)
{

	/* This subroutine sets up the initial conditions for Sod'd shock tube 
	 * problem assuming a constant value of the smoothing length.
	 */
	double pl	= 1.0;
	double denl	= 1.0;

	// conditions used by Herquist + Katz
	double denr	= 0.25;
	double pr	= 0.2154;

	// Zalezak problem
	//double denr	= 0.1;
	//double pr	= 0.01;

	for(int i=0; i<n; ++i) { partno[i] = i;}

	double delx		= (xmax-xmin)/(n-100);
	double delxl	= delx;
	double delxr	= delx;

	// for constant mass particles
	double mtot		= 0.75;
	
	int nleft=n/2;

	for(i=0; i<n; ++i) {
		
		if (i<=nleft) {
			
			delx	= delxl;
			den[i]	= denl;
			mass[i]	= den[i]*delx;
			vx[i]	= 0.0;
			if (i>1) {
				xx[i] = xx[i-1]+delx;
			}
			else {
				xx[i] = -delx*50.0;
			}
			hh[i]	= hfac*delx;
			p[i]	= pl;
			u[i]	= p[i]/den[i];
			u[i]	= u[i]/(gamma-1.0);
		}
		else {

			delx	= delxr;
			den[i]	= denr;
			mass[i]	= den[i]*delx;
			vx[i]	= 0.0;
			xx[i]	= xx[i-1]+delx;
			hh[i]	= hfac*delx;
			p[i]	= pr;
			u[i]	= p[i]/den[i];
			u[i]	= u[i]/(gamma-1.0);
		}

		if (ientro==1) u[i] = p[i]/pow(den[i],gamma);
	}

	xx[51] = 0.0;

}