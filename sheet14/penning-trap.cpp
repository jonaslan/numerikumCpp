#include <iostream>
#include "helfer.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "stepperbs.h"
#include <fstream>

using namespace std;

// Class implementing the differential equation (two-particle case) with the corresponding initial conditions
class Ode
{
	double B0;
	double freq2;
	double Omega;
	double Q2;
public:
	VecDoub ystart;
	Ode() : B0(1.), freq2(0.1), Omega(1.), Q2(0.1)
	{
		ystart.resize(8);
		// initial position
		ystart[0] = 0.3; // x1
		ystart[1] = 0.; // y1
		ystart[2] = 0.; // x2
		ystart[3] = -0.2; // y2

		// initial velocity
		ystart[4] = 1; //x1_dot
		ystart[5] = -1.; // y1_dot
		ystart[6] = 0.; // x2_dot
		ystart[7] = 0.; // y2_dot
	}

	void operator() (const double t, const VecDoub& y, VecDoub& dydt)
	{
		double r = sqrt((y[0]-y[2])*(y[0]-y[2])+(y[1]-y[3])*(y[1]-y[3]));
		dydt[0] = y[4];
		dydt[1] = y[5];
		dydt[2] = y[6];
		dydt[3] = y[7];
		dydt[4] = freq2*y[0]-Omega*y[5]+Q2*(y[0]-y[2])/r;
		dydt[5] = freq2*y[1]+Omega*y[4]+Q2*(y[1]-y[3])/r;
		dydt[6] = freq2*y[2]-Omega*y[7]-Q2*(y[0]-y[2])/r;
		dydt[7] = freq2*y[3]+Omega*y[6]-Q2*(y[1]-y[3])/r;
	}
};


int main ()
{
	ofstream file;
        file.open("penning-trap.txt",ios::trunc);

	// Instantiation of differential equation object
	Ode L;
	// Instantiation of two output objects, one for each solver
	Output out1(1000);
	//Output out2(1000);
	// Integration parameters
	const double tmax = 1e7;
	double atol = 1e-10;
	double rtol = 1e-10;
	double h1 = 0.01;
	double hmin = 0.;

	// Instantiation of two solver objects, one for each method 
	Odeint<StepperDopr853<Ode> > DoPrSolver(L.ystart, 0., tmax, atol, rtol, h1, hmin, out1, L);
	// Odeint<StepperBS<Ode> > BSSolver(L.ystart, 0., tmax, atol, rtol, h1, hmin, out2, L);

	// Solving the differential equation
	clock_t startDopr= clock();
	DoPrSolver.integrate();
	clock_t endDopr = clock();

	//clock_t startBS = clock();
	//	BSSolver.integrate();
	//clock_t endBS = clock();

	// Calculate #steps taken by each method: we didn't know if we should sum rejected and accepted steps, but that's what we did. We haven't observed differences anyway for the values used in the simulations.
	int steps_Dopr853 = out1.steps[0] + out1.steps[1];
	//int steps_BS = out2.steps[0] + out2.steps[1];
	
	// Print out number of steps and execution time taken by each method
	cout << "Dopr853: " << steps_Dopr853 << "\ttime: " << (double)(endDopr-startDopr)/(double)CLOCKS_PER_SEC << endl;
	//cout << "BS: " << steps_BS << "\ttime:" << (double)(endBS-startBS)/(double)CLOCKS_PER_SEC << endl;

	// Print out results stored in xsave (VecDoub) and ysave (MatDoub)
	for (int i=0; i < out1.count; ++i)
	{
		file << out1.xsave[i] << "\t" << out1.ysave[0][i] 
				      << "\t" << out1.ysave[1][i] 
				      << "\t" << out1.ysave[2][i] 
				      << "\t" << out1.ysave[3][i] 
				      << "\t" << out1.ysave[4][i] 
				      << "\t" << out1.ysave[5][i] 
				      << "\t" << out1.ysave[6][i] 
				      << "\t" << out1.ysave[7][i] << endl;
	}
	
	file.close();
}
