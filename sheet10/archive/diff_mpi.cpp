// TODO: Include required libraries, header files, and namespaces


// TODO: Integrand class goes here!


// **********
int main(int argc, char *argv[])
{
// 	Some necessary MPI variables
	int rank, size, worker_rank, idx=0;

// 	TODO: Initialize MPI environment


// 	TODO: Declare any additional variables


// 	Main branch: master process
	if (!rank)
	{
		int command, cnt=0;

		// TODO: Initialize output file. Find out where it needs to be closed!


		// Initialize master communication object
		Master m(size-1, MPI_DOUBLE, 2);
		while ( idx<N || m.some_working() )
		{
			// If a worker reports finished work, take results
			if (worker_rank = m.listen(&command))
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
							m.suspend_worker(worker_rank);
						break;

					case MW_return_result:
						m.get_result(worker_rank, ans);

						// TODO: Report reception of result and write to file


						break;

					case MW_job_done:
						m.free_worker(worker_rank);
						break;
				}
			usleep(100);
		}
	}

// 	Main branch: worker processes
	else
	{
		// TODO: Initialization of integration classes


		// Initialize worker communication object
		Worker w(MPI_DOUBLE, 2);

		// TODO: Repeatedly ask for work, do integration, and submit result to master


	}

//	Finalize the MPI environment
	MPI_Finalize();
}
