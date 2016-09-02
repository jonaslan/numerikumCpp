#include <iostream>
#include "helfer.h"
#include "odeint.h"
#include "stepperdopr853.h"
#include "stepperbs.h"
#include <omp.h>

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
	Lorentz L1;
	Lorentz L2;
	// Instantiation of two output objects, one for each solver
	Output out1(1000);
	Output out2(1000);
	// Integration parameters
	const double tmax = 1e7;
	const double atol = 1e-10;
	const double rtol = 1e-10;
	const double h1 = 0.01;
	const double hmin = 0.;

	// Specify number of threads
	omp_set_num_threads(2);

	// variables to measure execution time
	double startDopr, endDopr, startBS, endBS;
	double start = omp_get_wtime();
	// Parallelization
#pragma omp parallel sections
	{
#pragma omp section
		{
			Odeint<StepperDopr853<Lorentz> > DoPrSolver(L1.ystart, 0., tmax, atol, rtol, h1, hmin, out1, L1);
			startDopr= omp_get_wtime();
			DoPrSolver.integrate();
			endDopr = omp_get_wtime();
		}
#pragma omp section
		{
			Odeint<StepperBS<Lorentz> > BSSolver(L2.ystart, 0., tmax, atol, rtol, h1, hmin, out2, L2);
			startBS = omp_get_wtime();
			BSSolver.integrate();
			endBS = omp_get_wtime();
		}
	}
	double  end = omp_get_wtime();
	// Calculate #steps taken by each method: we didn't know if we should sum rejected and accepted steps, but that's what we did. We haven't observed differences anyway for the values used in the simulations.
	int steps_Dopr853 = out1.steps[0] + out1.steps[1];
	int steps_BS = out2.steps[0] + out2.steps[1];
	// Print out number of steps taken by each method
	cout << "Dopr853: " << steps_Dopr853 << "\ttime: " << (endDopr-startDopr) << endl;
	cout << "BS: " << steps_BS << "\ttime:" << (endBS-startBS) << endl;
	cout << (end - start) << endl;
}
