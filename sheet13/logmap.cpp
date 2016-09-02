#include <omp.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>

// Lyapunov exponent function
double lyapunov(double c, double* x)
{
	double sum = 0.;
	for (int i = 0; i < 10000; ++i)
	{
		sum += log(std::abs(c*(1.-2*x[i])));
	}
	return sum/10000.;
}

int main ()
{
	std::ofstream file1, file2;
	file1.open("chaos.txt", std::ios::trunc);
	file2.open("lyapunov.txt", std::ios::trunc);	

	// preparing parallel region
	int n = omp_get_num_procs();
	omp_set_num_threads(n);

	// coefficients array
	double c[10000];

	// splitting coefficient interval in 10000 parts, beginning with 1 and ending with 4
	for (int i = 0; i < 10000; ++i)
	{
		c[i] = 1 + 3./9999.*i;
	}

	double** A = new double*[10000];
	double* L = new double[10000];

	// Parallel region
	#pragma omp parallel shared(A,L)
	{
		// Parallel for loop
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < 10000; ++i)
		{
			// Thread-specific array
			double* values = new double[10000];

			// recursive for loop

			// here x should be declared 
			double x = 0.1;		
			for (int j = 0; j < 10000; ++j)
			{
				// and here's the really strange part. This implicit declaration of x works for some reason,
				// and stranger still is that without the implicit declaration the lyapunov exponent does not 
				// reach 1.4 or ln(2). Instead, if you don't have the explicit declaration, you get a bunch of
				// zeros in the mix which violates x(t) > 0 for all t. 
				x = c[i]*x*(1.-x);
				values[j] = x;
			}
			// exporting values
			A[i] = values;
		}
		
		// another parallel for loop for calculating the lyanpunov coefficients. Again, without an explicit declaration of
		// x above, the lyapunov exponent is 1.39 ~ 1.4 or ln(2) for c = 4. With the explicit declaration it ends up around 0.7
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < 10000; ++i)
		{
			double y = lyapunov(c[i],A[i]);
			L[i] = y;
		}
		
	}

	// printing the last hundred values values
	for (int i = 0; i < 10000; ++i)
	{
		double* t = A[i];
		for (int j = 9900; j < 10000; ++j)
			file1 << c[i] << "\t" << t[j] << std::endl;
		file2 << c[i] << "\t" << L[i] << std::endl;
	}


	// Here the arrays declared on the heap would be deleted if necessary. In this case the memory is freed on program exit
	file1.close();
	file2.close();	
}
