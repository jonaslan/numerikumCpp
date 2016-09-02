#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "adapt.h"
#include "masterworker.h"



// TODO: Integrand class goes here!
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
	int rank, size, worker_rank, idx=0;

// 	TODO: Initialize MPI environment

	MPI_Init(&argc,&argv);
	
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Status status;
	

// 	TODO: Declare any additional variables
	
	// integration array
	double ans[2];
	
	// integration resolution
	int N = 100;
	

// 	Main branch: master process
	if (!rank)
	{
		int command, cnt=0;

		// TODO: Initialize output file. Find out where it needs to be closed!

		ofstream file;
		file.open("mpi-integration.txt",std::ios::trunc);

		// Initialize master communication object
		Master m(size-1, MPI_DOUBLE, 2);
		
		while ( m.some_working() || idx < N)
		{
			if (worker_rank = m.listen(&command))
			{

				switch (command)
				{
					case MW_ask_for_job:
						if (idx < N)
						{
							m.send_work(worker_rank,&idx);
							++idx;
						}
						else
						{
							m.suspend_worker(worker_rank);
						}
						break;
					case MW_job_done:
						m.free_worker(worker_rank);
						break;
					case MW_return_result:
						
						m.get_result(worker_rank,&ans);
						file << std::setprecision(15) << ans[0] << "\t" << ans[1] << std::endl;
						break;
				}
			}
			usleep(100);
		}
		file.close();
	}

// 	Main branch: worker processes
	else
	{
		// TODO: Initialization of integration classes
		Diffusion D;
		Adapt A(1e-15);


		// Initialize worker communication object
		Worker w(MPI_DOUBLE, 2);

		// TODO: Repeatedly ask for work, do integration, and submit result to master
		int time;
		while( w.get_work( &time ) )
		{
			ans[0] = (double)time;
			D.set(ans[0]);
			ans[1] = A.integrate(D,1e-20,1.0);
			double integral = ans[1];
			w.send_result(&ans);
			w.done();
		}
	}
	MPI_Finalize();
}
