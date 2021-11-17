void sod(double xx[],double vx[],double mass[],double den[],
		 double p[],double u[],double hh[],double gamma,
		 double xmax,double xmin,int partno[],double hfac,
		 const int n,bool ientro)
{

	/* This subroutine sets up the initial conditions for Sod'd double
	 * shock tube problem.
	 */

    //changed//double* delx = new double[n];
	double delx;
      
	double pl	= 1.0;
    double denl = 1.0;

	// conditions used by Herquist + Katz
    double denr = 0.25;
    double pr	= 0.2154;

	// Zalezak problem
	//denr  = 0.1;
	//pr	= 0.01;

    for(int i=0; i<n; ++i) { partno[i]=i;}

	// for constant mass particles
	double	mtot = 0.75;
	int		neff = n-100;
	for(i=0; i<n; ++i) { mass[i] = mtot/neff;}

	double	delxl = mass[0]/denl;
	double	delxr = mass[0]/denr;

	double	xmid  = neff-(xmin/delxl)-(xmax/delxr);
			xmid  = xmid/((1.0/delxl)-(1.0/delxr));

	double	nleft = (1.0/delxl)*(xmid-xmin);
			nleft = nleft+50;

    for(i=0; i<n; ++i) { 
		if (i<=nleft) {
			delx	= delxl;
            den[i]	= denl;
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
            vx[i]	= 0.0;
            xx[i]	= xx[i-1]+delx;
            hh[i]	= hfac*delx;
            p[i]	= pr;
            u[i]	= p[i]/den[i];
            u[i]	= u[i]/(gamma-1.0);
		}

		// for constant size particles
		//hh(i)=hfac*delxr

		if (ientro==1) { u[i] = p[i]/pow(den[i],gamma);}

	}

    xx[51] = 0.0;
}
