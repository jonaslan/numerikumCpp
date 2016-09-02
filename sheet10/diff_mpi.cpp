// TODO: Include required libraries, header files, and namespaces
#include<mpi.h>
#include "masterworker.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include "adapt.h"
#include<ctime>

using namespace std;

// TODO: Integrand class goes here! (copied from diffusion.cpp provided by professor)

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
int main(int argc, char *argv[])
{
// 	Some necessary MPI variables
	int rank, size, worker_rank, idx=0; // idx will represent the current iteration 

// 	TODO: Initialize MPI environment
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//	Start timer
	clock_t start = clock();


// 	TODO: Declare any additional variables
	double ans[2];
	const int N = 200; // number of iterations
	const double tmax = 100.; // upper limit of integration
//      Determine time increment and smallest representable number
	const double dt = tmax/(N-1);
	const double eps = numeric_limits<double>::epsilon();

// 	Main branch: master process
	if (!rank)
	{
		int command, cnt=0;

// 	TODO: Initialize output file. Find out where it needs to be closed!
		ofstream data("diffusion_mpi.txt", ios::trunc);


// 	Initialize master communication object
		Master m(size-1, MPI_DOUBLE, 2); // argument with "2" is the max length of the array of data they pass to each other(here 2: time t and value of the integration, or idx)
		while ( idx<N || m.some_working() )
		{
//	If a worker reports finished work, take results
			if (worker_rank = m.listen(&command)) // checks m.listen() call for true (>0) or false (==0) and then assign worker_rank to the one listening
				switch (command)
				{
					case MW_ask_for_job:
						if (idx<N)
						{
							// Assign new task to worker
							m.send_work(worker_rank, &idx);
							++idx; 
						}
						else
							m.suspend_worker(worker_rank); /*this only applies when we have achieved the required name of steps but there is still the condition of the while which involves an OR and m.some_working() giving true. So we need to shut down all threads to get out of the loop. Could something similar have been achieved by a break sentence?*/
/* After this, the worker thread that had called w.get_work() and had sent a note of his waiting position (and was waiting for response), gets no response at all or empty response???? The difference  would be finishing the thread inside get_work or getting out of get_work!! What I think is that the main thread finishes (after shutting down or suspending every active worker) while the worker is still waiting, and when that happens an empty response is sent from the master to the worker, which then can get out of the get_work() and the associated loop and terminate correctly. Would that be correct? */
						break;

					case MW_return_result:
						m.get_result(worker_rank, ans);

						// TODO: Report reception of result and write to file. 

						++cnt; // Is this the report of reception? Other places to use it? Write to the file??
						data << ans[0] << '\t' << ans[1] << endl; //write t and integral estimate to file 


						break;

					case MW_job_done:
						m.free_worker(worker_rank); //once job is done, the master must update its array containing the states of each worker and declare it as free
						break;
				}
			usleep(100);
		}

		// Close file
		data.close();
	}

// 	Main branch: worker processes
	else
	{
// TODO: Initialization of integration classes --> Only one!! Maybe ans should be defined here as well?
		Adapt adapt(1e-15);
		Diffusion f;


// Initialize worker communication object
		Worker w(MPI_DOUBLE, 2);

// TODO: Repeatedly ask for work, do integration, and submit result to master

		while (w.get_work(&idx))
		{

			double t = idx * dt; // calculate time corresponding to iteration
			f.set(t); // set time in the integrand
			ans[0] = t;
			ans[1] = adapt.integrate(f, eps, 1.); // Remember: 1. is selected as upper bound because the integral is improper and we had to do the "trick"
			w.send_result(ans); //send result back to the master

			w.done(); // It must wait until the master is "willing" to listen its message. When master does so, master updates his state of working workers
		}
	}
// 	Finish timer
clock_t end = clock();

// Show a list in standard output with the execution time of each process
cout << "CPU time (in sec) taken by process " << rank << ": " << (end - start)/(double)CLOCKS_PER_SEC << endl;

//	Finalize the MPI environment
	MPI_Finalize();
}
