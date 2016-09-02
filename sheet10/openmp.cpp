#include<omp.h>
#include<iostream>
#include<unistd.h>
#include<ctime>
#include<fstream>


int main ()
{
	int sum = 0; // If we define it here, the variable is shared in every loop (unless we specify the contrary, which is not the case)
	// int sum[2] = {0,0}; // This is to check how the result of the for loop is different for different threads if we use a different summation variable (sum[0],sum[1]) for each thread
	
	
	// Open file that stores which iterations have been assigned to each thread 
	std::ofstream file;
	file.open("parallel_test.txt", std::ios::trunc);
	
	// Conditional compilation directives allow us to easily switch from serial to parallel execution
	#ifdef _OPENMP
	omp_set_num_threads( omp_get_num_procs() ); // Set total number of threads from the total number of processors in the machine
	int nt; // To store number of threads

	#pragma omp parallel shared(sum) // Start of parallel region
	#endif
	{
		#ifdef _OPENMP		
		nt = omp_get_num_threads(); // It only makes sense to use it (to get a result different from 1) inside the parallel region. 
		#pragma omp for schedule(dynamic) // distributed for loop with a dynamic scheduling
		#endif
		for (int i = 0; i < 100; ++i)
		{
			#ifdef _OPENMP
			#pragma omp critical // This ensures only one thread at a time attempts to write to the file
			
			file << i << "\t" << "executed in thread #" << omp_get_thread_num() << std::endl;
			#endif
	
			sum += i; // We don't run into problems due to the shared nature of the variable since addition is commutative
			//sum[omp_get_thread_num()] += i;
		}
	}
	
	std::cout << "Sum of [0,99] = " << sum  << " calculated on " << nt << " processors" << std::endl;
	std::cout << "Sum of [0,99]  calculated through formula of sum of finite arithmetic series: 100*(99-0)/2 = 4950" << std::endl;
	// std::cout << sum[0] << " " << sum[1] << std::endl; // shows that the two variables are in general different (they depend on the thread)
	
	file.close();
}
