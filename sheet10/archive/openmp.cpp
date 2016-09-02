#include<omp.h>
#include<iostream>
#include<unistd.h>
#include<ctime>
#include<fstream>



int main ()
{
	//int sum = 0;
	int sum[2] = {0,0};
	
	
	
	std::ofstream file;
	file.open("Parallel test.txt", std::ios::trunc);
	
	#ifdef _OPENMP
	omp_set_num_threads( omp_get_num_procs() );
	int nt = omp_get_num_threads();
		
	#pragma omp parallel shared(sum)
	#endif
	{
		#ifdef _OPENMP		
		#pragma omp for schedule(dynamic)
		#endif
		for (int i = 0; i < 100; ++i)
		{
			#ifdef _OPENMP
			#pragma omp critical
			
			//file << i << "\t" << omp_get_thread_num() << std::endl;
			#endif
	
			sum[omp_get_thread_num()] += i;
		}
	}
	
	//std::cout << "Sum of [0,99] = " << sum  << " calculated on " << nt << " processors" << std::endl;
	
	std::cout << sum[0] << " " << sum[1] << std::endl;
	
	file.close();
}