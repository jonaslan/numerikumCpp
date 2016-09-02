#include "helfer.h"
#include "adapt.h"
#include <limits>
#include <iostream>
#include <fstream>
#include<omp.h>
#include <ctime>

// ******* Integrand class (copied from professor's solution of previous exercise sheet)
class Diffusion
{
	double lambda, border;
	double s, t;

// 	Inline auxiliary function for the spectrum g
	inline double g(const double x)
	{
		return pow( 1.+x*x, -s/2. );
	}

// 	Auxiliary function for the integrand
	double integrand(const double x)
	{
		// Define constants
		double f1, f2, f3;
		const double f4 = exp( -t/(2.*lambda) );

		// Select range
		if (x < border)
		{
			f1 = sqrt( 1./(4.*lambda*lambda) - x*x/3. );
			f2 = 1./f1;
			f3 = sinh(t*f1);
			return g(x) * f2 * f3 * f4;
		}
		else if (x > border)
		{
			f1 = sqrt( -1./(4.*lambda*lambda) + x*x/3. );
			f2 = 1./f1;
			f3 = sin(t*f1);
			return g(x) * f2 * f3 * f4;
		}
		else
			return g(x) * t * f4;
	}

	public:
// 	Constructor: initialize parameters
	Diffusion()
	{
		lambda = 3.;
		border = sqrt(3.)/(2.*lambda);
		s = 5./3.;
		t = 0.;
	}

// 	Set the time for each integration
	void set(const double tt)
	{
		t = tt;
	}

// 	Return combined integrand for x<1 and for x>1
	double operator() (const double x)
	{
		const double y = 1./x;
		return integrand(x) + y*y * integrand(y);
	}
};

// **********
int main()
{
//	Set number of time steps and maximum time
	const int N = 200;
	const double tmax = 100.;

// 	Determine time increment and smallest representable number
	const double dt = tmax/(N-1);
	const double eps = numeric_limits<double>::epsilon();


// 	Instantiate integrand object
	Diffusion f;

// 	Instantiate adaptive integration object
	Adapt adapt(1e-15);

//	Declare variable that store the time t at which integration is performed and the corresponding numerical approximation
	double t, ans;

//	Set number of threads depending on available processors and compilation flag _OPENMP. Files to store output are opened	
#ifdef _OPENMP
	omp_set_num_threads( omp_get_num_procs() );
	cout << omp_get_num_procs() << endl;
	ofstream datei("diffusion_omp.txt",ios::trunc);
#else
	ofstream datei("diffusion_non-parallel.txt",ios::trunc);
#endif
	
#ifdef _OPENMP
//	Timer start for parallel execution (gives result directly in seconds). It is not possible to calculate total execution time through the normal clock() function since this gives total CPU time, i.e. the sum of CPU time of each thread
	double start = omp_get_wtime();
#pragma omp parallel
{
//	Start distributed for with the variables being used inside each iteration defined as private
#pragma omp for private(t, ans, f, adapt) //schedule(dynamic) //possible option: it doesn't change much the running time in our machine
#else
	clock_t start = clock(); // Start timer in case of non-parallel execution (results in CPU clocks)
#endif
	for (int i=0; i<N; ++i)
	{
		// Set time in integrand
		t = i*dt;
		f.set(t);

		// Perform integration for time t
		ans = adapt.integrate(f,eps,1.);

		// Output result as a function of time
#ifdef _OPENMP
#pragma omp critical // Avoid multiple threads trying to write to a file at the same time
#endif
		datei << t << '\t' << ans << endl; // writing to file

		// Warn if integration fails
		if (adapt.offTolerance())
			cout << "Required tolerance has not been met" << endl;
	}

#ifdef _OPENMP
}
double end = omp_get_wtime(); //End time for parallel execution
#else
clock_t end = clock(); //End time (in clocks) for serial execution
#endif
	datei.close(); //close output file

	// Calculate execution time in seconds
	std::cout << (end-start)
	#ifndef _OPENMP
	/(double)CLOCKS_PER_SEC
	#endif
	<< std::endl;
}
