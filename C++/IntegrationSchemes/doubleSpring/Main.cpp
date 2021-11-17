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

	double f   = 10.0;	
	double x01 = 0.0;
	double v01 = 0.0;
	double x02 = 0.00;
	double v02 = 0.0;
	double t = 0.0;
	double dt = 1e-4;

	double m1 = 1.0;
	double k1 = 10000.0;
	double c1 = 200.0;
	double x1 = x01;
	double v1 = v01;
	double a1 = 0.0;
	double m2 = 1.0;
	double k2 = 10000.0;
	double c2 = 200.0;
	double x2 = x02;
	double v2 = v02;
	double a2 = 0.0;
	double LF_x1_t = x01;
	double LF_x1_tpdt = x01;
	double LF_v1_t = v01;
	double LF_v1_tmdt = v01;
	double LF_v1_tpdt2 = v01;
	double LF_v1_tmdt2 = v01;
	double LF_a1_t = 0.0;
	double LF_a1_tmdt = 0.0;
	double LF_x2_t = x02;
	double LF_x2_tpdt = x02;
	double LF_v2_t = v02;
	double LF_v2_tmdt = v02;
	double LF_v2_tpdt2 = v02;
	double LF_v2_tmdt2 = v02;
	double LF_a2_t = 0.0;
	double LF_a2_tmdt = 0.0;

	// exact solutions
	double omega_n1 = sqrt(k1/m1);
	double omega_n2 = sqrt(k2/m2);
	double ksi1 = c1/(2.0*sqrt(m1*k1));
	double ksi2 = c2/(2.0*sqrt(m2*k2));
	double omega_d1 = omega_n1*sqrt(1.0-(ksi1*ksi1));
	double omega_d2 = omega_n2*sqrt(1.0-(ksi2*ksi2));

	cout<<"ksi1: "<<ksi1<<endl;
	cout<<"ksi2: "<<ksi2<<endl;
	cout<<"omega_n1: "<<omega_n1<<endl;
	cout<<"omega_n2: "<<omega_n2<<endl;
	cout<<"omega_d1: "<<omega_d1<<endl;
	cout<<"omega_d2: "<<omega_d2<<endl;
	double tend = 1;
	for(t=0; t<=tend; )
		//for(int i=0; i<100001; ++i)
	{
//		double amp = 10.0; // amplitude (m)
//		double freq = 1.0; // frequency (rad/s)
//		f = amp*sin(freq*t);
		// new acceleration
		LF_v2_t = 0.25*(LF_v2_tmdt+dt*LF_a2_tmdt) + 0.75*(LF_v2_tmdt2+0.5*dt*LF_a2_tmdt);
		LF_v1_t = 0.25*(LF_v1_tmdt+dt*LF_a1_tmdt) + 0.75*(LF_v1_tmdt2+0.5*dt*LF_a1_tmdt); // extrapolation for v(t)
		LF_a1_t = (-k1*LF_x1_t-c1*LF_v1_t+k2*(LF_x2_t-LF_x1_t)+c2*(LF_v2_t-LF_v1_t))/m1; // a(t)=f(x(t),v(t))
		LF_a2_t = (-k2*(LF_x2_t-LF_x1_t)-c2*(LF_v2_t-LF_v1_t)+f)/m2;
		a1 = (-k1*x1-c1*v1+k2*(x2-x1)+c2*(v2-v1))/m1;
		a2 = (-k2*(x2-x1)-c2*(v2-v1)+f)/m2;

		// move
		LF_v1_tpdt2 = LF_v1_tmdt2 + dt*LF_a1_t; // v(t+dt/2)=v(t-dt/2)+dt*a(t)
		LF_x1_tpdt = LF_x1_t + LF_v1_tpdt2*dt; // x(t+dt)=x(t)+dt*v(t+dt/2)
		LF_v2_tpdt2 = LF_v2_tmdt2 + dt*LF_a2_t;
		LF_x2_tpdt = LF_x2_t + LF_v2_tpdt2*dt;
		v1 = v1 + a1*dt;
		v2 = v2 + a2*dt;
		x1 = x1 + v1*dt;
		x2 = x2 + v2*dt;


		// save
		outFile<<t
				  <<" "<<x1<<" "<<v1<<" "<<a1
			      <<" "<<x2<<" "<<v2<<" "<<a2
				  <<" "<<LF_x1_t<<" "<<LF_v1_t<<" "<<LF_a1_t
			      <<" "<<LF_x2_t<<" "<<LF_v2_t<<" "<<LF_a2_t
				<<endl;
		
		// reset
		a1 = 0.0; 
		a2 = 0.0;
		
		//proceed
		t = t + dt;
		LF_x1_t = LF_x1_tpdt;			// x(t+dt)	->x(t)
		LF_v1_tmdt2 = LF_v1_tpdt2;	// v(t+dt/2)->v(t-dt/2)
		LF_v1_tmdt = LF_v1_t;			// v(t)		->v(t-dt)
		LF_a1_tmdt = LF_a1_t;			// a(t)		->a(t-dt)
		LF_x2_t = LF_x2_tpdt;			
		LF_v2_tmdt2 = LF_v2_tpdt2;	
		LF_v2_tmdt = LF_v2_t;		
		LF_a2_tmdt = LF_a2_t;		
	}

	outFile.close();
	return 0;
}