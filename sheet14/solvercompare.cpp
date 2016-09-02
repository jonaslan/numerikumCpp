#include <iostream>
#include "helfer.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "stepperbs.h"

using namespace std;

// Class implementing the differential equation with the corresponding initial conditions
class Lorentz
{
	double B0;
public:
	VecDoub ystart;
	Lorentz() : B0(1.)
	{
		ystart.resize(6);
		// initial position
		ystart[0] = 0.;
		ystart[1] = 0.;
		ystart[2] = 0.;

		// initial velocity
		ystart[3] = 1./sqrt(2);
		ystart[4] = 0;
		ystart[5] = 1./sqrt(2);
	}

	void operator() (const double t, const VecDoub& y, VecDoub& dydt)
	{
		dydt[0] = y[3];
		dydt[1] = y[4];
		dydt[2] = y[5];
		dydt[3] = y[4]*B0;
		dydt[4] = -y[3]*B0;
		dydt[5] = 0.;
	}
};


int main ()
{
	// Instantiation of differential equation object
	Lorentz L;
	// Instantiation of two output objects, one for each solver
	Output out1(1000);
	Output out2(1000);
	// Integration parameters
	const double tmax = 1e7;
	double atol = 1e-10;
	double rtol = 1e-10;
	double h1 = 0.01;
	double hmin = 0.;
	// Instantiation of two solver objects, one for each method 
	Odeint<StepperDopr853<Lorentz> > DoPrSolver(L.ystart, 0., tmax, atol, rtol, h1, hmin, out1, L);
	Odeint<StepperBS<Lorentz> > BSSolver(L.ystart, 0., tmax, atol, rtol, h1, hmin, out2, L);
	// Solving the differential equation
clock_t startDopr= clock();
	DoPrSolver.integrate();
clock_t endDopr = clock();
clock_t startBS = clock();
	BSSolver.integrate();
clock_t endBS = clock();
	// Calculate #steps taken by each method: we didn't know if we should sum rejected and accepted steps, but that's what we did. We haven't observed differences anyway for the values used in the simulations.
	int steps_Dopr853 = out1.steps[0] + out1.steps[1];
	int steps_BS = out2.steps[0] + out2.steps[1];
	// Print out number of steps taken by each method
	cout << "Dopr853: " << steps_Dopr853 << "\ttime: " << (double)(endDopr-startDopr)/(double)CLOCKS_PER_SEC << endl;
	cout << "BS: " << steps_BS << "\ttime:" << (double)(endBS-startBS)/(double)CLOCKS_PER_SEC << endl;

}
