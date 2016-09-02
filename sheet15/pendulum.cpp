#include <iostream>
#include <cmath>
#include <ctime>
#include "helfer.h"
#include "stepperdopr853.h"
#include "stepperstoerm.h"
#include "odeint.h"


using namespace std;

class Pendulum
{
	double g;
	double r;

public:
	
	VecDoub ystart;	
	
	Pendulum(const double rIn, const double theta0In) : g(9.81), r(rIn) 
	{
		ystart.resize(2);
		// initial position
		ystart[0] = theta0In;
		
		// initial velocity
		ystart[1] = 0;
	}
	
	void operator () (const double t, const VecDoub& z, VecDoub& dydt)
	{
		dydt[0] = -g/r*sin(y[0]);
		dydt[1] = y[1];
	}
};


int main ()
{

	double r = 1.;
	double theta0 = 0.5*M_PI;
	
	Pendulum P(r,theta0);
	Output out1(1000);
	Output out2(1000);
	// Integration parameters
	const double tmax = 100.;
	double atol = 1e-10;
	double rtol = 1e-10;
	double h1 = 0.01;
	double hmin = 0.;
	// Instantiation of two solver objects, one for each method 
	Odeint<StepperDopr853<Pendulum> > DoPrSolver(P.ystart, 0., tmax, atol, rtol, h1, hmin, out1, P);
	Odeint<Stepperstoerm<Pendulum> > StoermSolver(P.ystart, 0., tmax, atol, rtol, h1, hmin, out2, P);
	// Solving the differential equation
	clock_t dprstart = clock();
	DoPrSolver.integrate();
	clock_t	dprend = clock();
	
	clock_t stoermstart = clock();
	StoermSolver.integrate();
	clock_t stoermend = clock();
	
	
	// Calculate #steps taken by each method: we didn't know if we should sum rejected and accepted steps, but that's what we did. We haven't observed differences anyway for the values used in the simulations.
	int steps_Dopr853 = out1.steps[0] + out1.steps[1];
	int steps_Stoerm = out2.steps[0] + out2.steps[1];
	
	cout << "Dormand-Prince: " << dprend-dprstart << " steps: "steps_Dopr853 << "\n " << steps_stoerm << endl;
	
}