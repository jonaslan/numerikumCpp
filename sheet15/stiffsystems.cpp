// stiff systems of differential equations
#include <iostream>
#include "helfer.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "stepperstoerm.h"
#include "stepperross.h"
#include "steppersie.h"


using namespace std;


class Robertson
{

public:
	VecDoub ystart;
	Robertson()
	{
		ystart.resize(3);
		// initial position
		ystart[0] = 1.;
		ystart[1] = 1.;
		ystart[2] = 0.;
	}
	
	void operator() (const double t, const VecDoub& y, VecDoub& dydt)
	{
		dydt[0] = -0.013*y[0]-1000.*y[0]*y[2];
		dydt[1] = -2500.*y[1]*y[2];
		dydt[2] = -0.013*y[0]-1000.*y[0]*y[2]-2500.*y[1]*y[2];
	}
	
	void jacobian (const double x, const VecDoub& y, VecDoub& dfdx, MatDoub& dfdy)
	{
		dfdx[0] = 0.;
		dfdx[1] = 0.;
		dfdx[2] = 0.;
		
		dfdy[0][0] = -0.013 - 1000.*y[2];
		dfdy[0][1] = 0.;
		dfdy[0][2] = -1000.*y[0];
		dfdy[1][0] = 0.;
		dfdy[1][1] = -2500.*y[2];
		dfdy[1][2] = -2500.*y[1];
		dfdy[2][0] = -0.013 - 1000.*y[2];
		dfdy[2][1] = -2500.*y[2];
		dfdy[2][2] = -1000.*y[0]-2500.*y[1];	
	}
};


int main ()
{

	Robertson R;
	Output out1(1000);
	Output out2(1000);
	Output out3(1000);
	Output out4(1000);
	
	// Integration parameters
	const double tmax = 1e5;
	double atol = 1e-6;
	double rtol = 1e-6;
	double h1 = 2.9e-4;
	double hmin = 0.;
	
	// Instantiation of two solver objects, one for each method 
	
	Odeint<StepperDopr853<Robertson> > DoPrSolver(R.ystart, 0., tmax, atol, rtol, h1, hmin, out1, R);
	Odeint<StepperStoerm<Robertson> > StoermSolver(R.ystart, 0., tmax, atol, rtol, h1, hmin, out2, R);
	
	Odeint<StepperRoss<Robertson> > RossSolver(R.ystart, 0., tmax, atol, rtol, h1, hmin, out3, R);
	Odeint<StepperSie<Robertson> > SieSolver(R.ystart, 0., tmax, atol, rtol, h1, hmin, out4, R);
	// Solving the differential equation
	clock_t dprstart = clock();
	DoPrSolver.integrate();
	clock_t	dprend = clock();
	
	clock_t stoermstart = clock();
	StoermSolver.integrate();
	clock_t stoermend = clock();
	
	clock_t rossstart = clock();
	RossSolver.integrate();
	clock_t rossend = clock();
	
	clock_t siestart = clock();
	SieSolver.integrate();
	clock_t sieend = clock();
	
	
	// Calculate #steps taken by each method: we didn't know if we should sum rejected and accepted steps, but that's what we did. We haven't observed differences anyway for the values used in the simulations.
	int steps_Dopr853 = out1.steps[0];
	int steps_Stoerm = out2.steps[0];
	int steps_Ross = out3.steps[0];
	int steps_Sie = out4.steps[0];
	cout << "Ross: " << rossend-rossstart << " steps: " << steps_Ross << "\n " << steps_Sie << endl;

	
// 	cout << "Dormand-Prince: " << dprend-dprstart << " steps: " << steps_Dopr853 << "\n " << steps_Stoerm << endl;
}