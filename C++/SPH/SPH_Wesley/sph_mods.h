/*********************************************************************
 * sph_mods.h header file for the sph_mods.cc subroutines. By Wesley
 * Miaw <wesley@wesman.net>. Original sph_mods.f code by Kevin Olson
 * <olson@jeans.gsfc.nasa.gov> available from NASA/Goddard at
 * http://sdcd.gsfc.nasa.gov/ESS/exchange/contrib/olson/sph_1d.html.
 *
 * Modified: 05/08/2001
 ********************************************************************/
#include <math.h>

void accelerations(float *xx, float *vx, float *ax, float *yy, float
		   *vy, float *ay, float *zz, float *vz, float *az,
		   float *au, float *hh, float *p, float *u, float
		   *mass, float *den, float *divv_x, float *divv_y,
		   float *divv_z, float *sound, float *maxmu_x, float
		   *maxmu_y, float *maxmu_z, int n, int inorm, float
		   alpha, float beta, int ibulk, float alphab, float
		   betab, float gamma, int ientro, int ientro2);

void bouyancy(float *xx, float *yy, float *zz, float *az, float *den,
	      int n, float b1, float grav);

void therm_diffusion(float *xx, float *vx, float *yy, float *vy, float
		     *zz, float *vz, float *hh, float *mass, float
		     *den, float *p, float *u, float *divv_x, float
		     *divv_y, float *divv_z, float *sound, int n,
		     float g1, float g2, float gamma, float
		     *thermdiff);

void push1(float *xx, float *yy, float *zz, float *u, float *vx, float
	   *vy, float *vz, float *ax, float *ay, float *az, float *au,
	   float c1, float c2, int n);

void benz_h1(float *hh, float *divv, float c1, int n);

void push2(float *xx, float *u, float *vx, float *ax, float *au, float
	   c1, float c2, float delt, float *xo, float *vxo, float *uo,
	   float t, int n);

void benz_h2(float *hh, float c1, float delt, float *divv, float *ho,
	     float t, int n);

void update_den(float *xx, float *yy, float *zz, float *mass, float
		*den, float *hh, int n);

void update_p(float *hh, float *mass, float *u, float *den, float *p,
	      float gamma, int n, int ientro);

void new_delt(float *hh, float *divv_x, float *divv_y, float *divv_z,
	      float *sound, int inorm, float alpha, float beta, int
	      ibulk, float alphab, float betab, float *maxmu_x, float
	      *maxmu_y, float *maxmu_z, int n, float c1, float c2,
	      float cour, float delt);

void reflecting(float *xx, float *vx, float *xo, float *vxo, float
		*yy, float *vy, float *yo, float *vyo, float *zz,
		float *vz, float *zo, float *vzo, float *den, float
		*mass, float *p, float *u, float *uo, float *hh, float
		*ho, int n);

void ww(float *r, float *hh, float *ht, int n, float *w);

void gradww(float *r, float *hh, float *ht, int n, float *gradw);

void smooth(float *xx, float *yy, float *zz, float *hh, float *mass,
	    float *den, int n, float *u);

void div_v(float *xx, float *vx, float *yy, float *vy, float *zz,
	   float *vz, float *hh, float *mass, float *den, int n,
	   float *divv_x, float *divv_y, float *divv_z);

void SM_h(float *xx, float *yy, float *zz, float *hh, float *mass,
	  float *den, float *divv_x, float *divv_y, float *divv_z,
	  float t, float delt, int n);

void smooth_den(float *xx, float *yy, float *zz, float *hh, float
		*mass, float *den, int n, float *densmooth);

void d_smooth_den(float *xx, float *yy, float *zz, float *hh, float
		  *mass, float *den, float *densmooth, float *divv_x,
		  float *divv_y, float *divv_z, int n, float
		  *dsmoothdt);

bool any_overlap(float *r, float *hh, float *ht, int n);

void lshift(float *a, int n);

void rshift(float *a, int n);

float minval(float *a, int n);

float abs(float x);

float min(float x, float y);
