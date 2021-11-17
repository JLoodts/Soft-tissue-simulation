#include <math.h>
#pragma warning(disable:4244)
#pragma warning(disable:4305)
void sod(float *xx, float *vx, float *yy, float *vy, float *zz,
	 float *vz, float *mass, float *den, float *p, float *u,
	 float *hh, float gamma, float xmax, float xmin, int
	 *partno, float hfac, int n, int ientro)
{
  
  float delx;
  float denl, denr, pl, pr, mtot;

  float delxl, delxr, xmid;
  int neff, nleft;

  pl = 1.0;
  denl = 1.0;

  // conditions used by herquist and katz
  
  denr = 0.25; 
  pr = 0.2154;

  for (int i=0; i<n; i++) 
    partno[i]=i+1;

  // for constant mass particles

  mtot = 0.75;
  neff = n-100; // ?
  
  for (i=0; i<n; i++)
    mass[i] = mtot/neff;  //?
  
  delxl = mass[1]/denl;
  delxr = mass[1]/denr;

  xmid = neff - (xmin/delxl)-(xmax/delxr);
  xmid = xmid/((1.0/delxl)-(1.0/delxr));

  nleft = (int)((1.0/delxl)*(xmid-xmin));
  nleft = nleft + 50;
  
  int l=0;

  for (i=0; i<n; i++) {
    if (i<=nleft) {
      delx = delxl;
      den[i]=denl;
      vx[i]=0.0;
      if (i>l)
	xx[i]=xx[i-1]+delx;
      else
	xx[i]=-delx*50.0;
      hh[i] = hfac * delx;
      p[i] = pl;
      u[i] = p[i]/den[i];
      u[i] = u[i]/(gamma-1.0);
    }
    yy[i] = 0; vy[i] = 0;
    zz[i] = 0; vz[i] = 0;
    
    // for constant size particles
    // hh[i] = hfac*delxr;
    
    if (ientro == 1)
      u[i]=p[i]/pow(den[i],gamma);
  }

  xx[50] = 0.0;

}
