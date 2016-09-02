#include <omp.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>

double lyapunov(double c, double* x)
{
	double sum = 0.;
	for (int i = 0; i < 10000; ++i)
	{
		sum += log(std::abs(c*(1.-x[i])));
	}
	return sum/10000.;
}

int main ()
{
	std::ofstream file1, file2;
	file1.open("chaos.txt", std::ios::trunc);
	file2.open("lyapunov.txt", std::ios::trunc);	

	int n = omp_get_num_procs();

	omp_set_num_threads(n);

	double c[10000];

	for (int i = 0; i < 10000; ++i)
	{
		c[i] = 1 + 3./9999.*i;
	}

	double** A = new double*[10000];
	double* L = new double[10000];

	#pragma omp parallel shared(A)
	{
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < 10000; ++i)
		{
			double* values = new double[10000];
			for (int j = 0; j < 10000; ++j)
			{
				double x = c[i]*x*(1.-x);
				values[j] = x;
			}
			A[i] = values;
		}
		
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < 10000; ++i)
		{
			double x = lyapunov(c[i],A[i]);
			L[i] = x;
		}
		
	}

	for (int i = 0; i < 10000; ++i)
	{
		double* t = A[i];
		for (int j = 9900; j < 10000; ++j)
			file1 << c[i] << "\t" << t[j] << std::endl;
		file2 << c[i] << "\t" << L[i] << std::endl;
	}


	file1.close();
	file2.close();	
}
