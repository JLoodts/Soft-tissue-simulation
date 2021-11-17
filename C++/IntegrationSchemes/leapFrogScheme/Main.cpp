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

	double x01 = 0.01;
	double v01 = 0.0;
	double x02 = 0.0;
	double v02 = 0.0;
	double t = 0.0;
	double dt = 1e-4;

	double m1 = 1.0;
	double k1 = 10000.0;
	double c1 = 200.0;
	double x1 = x01;
	double v1 = v01;
	double v1minus = v01;
	double v1plus = v01;
	double a1 = 0.0;
	double m2 = 1.0;
	double k2 = 10000.0;
	double c2 = 200.0;
	double x2 = x02;
	double v2 = v02;
	double v2minus = v02;
	double v2plus = v02;
	double a2 = 0.0;

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
	double tend = 0.2;
	for(t=0; t<=tend; )
		//for(int i=0; i<100001; ++i)
	{
		// new acceleration
		a1 = (-k1*x1-c1*v1+k2*(x2-x1)+c2*(v2-v1))/m1;
		a2 = (-k2*(x2-x1)-c2*(v2-v1))/m2;

		// move
		v1plus = v1minus + a1*dt;
		v2plus = v2minus + a2*dt;
		x1 = x1 + v1plus*dt;
		x2 = x2 + v2plus*dt;
		v1 = (v1plus + v1minus)*0.5;
		v2 = (v2plus + v2minus)*0.5;

		// save
		outFile<<t<<" "<<x1<<" "<<v1<<" "<<a1
			      <<" "<<x2<<" "<<v2<<" "<<a2<<endl;
		
		// reset
		a1 = 0.0; 
		a2 = 0.0;
		
		//proceed
		t = t + dt;
		v1minus = v1plus;
		v2minus = v2plus;
	}

	outFile.close();
	return 0;
}