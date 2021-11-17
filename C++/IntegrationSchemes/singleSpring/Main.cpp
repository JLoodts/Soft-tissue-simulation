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

// the main function
// -----------------
int main(int argc, char** argv)
{
	const char* nameOut = "data/forceDeformation.txt";
	ofstream outFile;
	outFile.open(nameOut, ios::out);

	double x0 = 0.01;
	double v0 = 0.0;
	double t = 0.0;
	double dt = 5e-4;

	double m = 1.0;
	double k = 1000.0;
	double c = 200.0;
	double x = x0;
	double v = v0;
	double a = 0.0;
	double LF_x_t = x0;
	double LF_x_tpdt = x0;
	double LF_v_t = v0;
	double LF_v_tmdt = v0;
	double LF_v_tpdt2 = v0;
	double LF_v_tmdt2 = v0;
	double LF_a_t = 0.0;
	double LF_a_tmdt = 0.0;

	// exact solutions
	double omega_n = sqrt(k/m);
	double ksi = c/(2.0*sqrt(m*k));
	double omega_d = omega_n*sqrt(1.0-(ksi*ksi));
	double x_e = x0;
	//double v_e = v0;
	//double a_e = 0.0;

	cout<<"ksi: "<<ksi<<endl;
	cout<<"omega_n: "<<omega_n<<endl;
	cout<<"omega_d: "<<omega_d<<endl;
/**/	double max_error = 0;
	double LF_max_error = 0;

	for(dt = 1e-7; dt<3e-3; dt *=1.01)
	{
		for(t=0; t<=0.1+dt; )
		//for(int i=0; i<100001; ++i)
		{
			// new acceleration
			a = (-k*x-c*v)/m;

			LF_v_t = 0.25*(LF_v_tmdt+dt*LF_a_tmdt) + 0.75*(LF_v_tmdt2+0.5*dt*LF_a_tmdt); // extrapolation for v(t)
			LF_a_t = (-k*LF_x_t-c*LF_v_t)/m; // a(t)=f(x(t),v(t))

			// move
			v = v + a*dt;
			x = x + v*dt;
			
			LF_v_tpdt2 = LF_v_tmdt2 + dt*LF_a_t; // v(t+dt/2)=v(t-dt/2)+dt*a(t)
			LF_x_tpdt = LF_x_t + LF_v_tpdt2*dt; // x(t+dt)=x(t)+dt*v(t+dt/2)
						

			if(ksi < 1.0) {
				x_e = exp(-ksi*omega_n*t)*(x0*cos(omega_d*t)+(v0+ksi*omega_n*x0)/omega_d*sin(omega_d*t));
			} else {
				if(ksi > 1.0) {
					double sqrtKsi2One = sqrt(ksi*ksi-1);
					x_e = (x0*omega_n*(ksi+sqrtKsi2One)+v0)/(2.0*omega_n*sqrtKsi2One)*exp((-ksi+sqrtKsi2One)*omega_n*t)+
						 (-x0*omega_n*(ksi-sqrtKsi2One)-v0)/(2.0*omega_n*sqrtKsi2One)*exp((-ksi-sqrtKsi2One)*omega_n*t);
				} else {
					x_e = exp(-omega_n*t)*(x0+(v0+omega_n*x0)*t);
				}
			}
			// save
			//outFile<<t<<" "<<x<<" "<<v<<" "<<a<<" "<<x_e
			//		  <<" "<<LF_x_t<<" "<<LF_v_t<<" "<<LF_a_t<<endl;
/**/			if(fabs(x-x_e)>max_error)max_error = x-x_e;
			if(fabs(LF_x_t-x_e)>LF_max_error)LF_max_error = LF_x_t-x_e;

			// reset
			a = 0.0;
			//a_e = 0.0;
			
			//proceed
			t = t + dt;
			LF_x_t = LF_x_tpdt;			// x(t+dt)	->x(t)
			LF_v_tmdt2 = LF_v_tpdt2;	// v(t+dt/2)->v(t-dt/2)
			LF_v_tmdt = LF_v_t;			// v(t)		->v(t-dt)
			LF_a_tmdt = LF_a_t;			// a(t)		->a(t-dt)
		}
		// save
/**/		outFile<<dt<<" "<<max_error<<" "<<LF_max_error<<" "<<endl;
		
		max_error = 0;
		LF_max_error = 0;
		
		t=0;
		x = x0;
		v = v0;

		LF_x_t = x0;
		LF_x_tpdt = x0;
		LF_v_t = v0;
		LF_v_tmdt = v0;
		LF_v_tpdt2 = v0;
		LF_v_tmdt2 = v0;
		LF_a_t = 0.0;
		LF_a_tmdt = 0.0;
	}

	outFile.close();
	return 0;
}