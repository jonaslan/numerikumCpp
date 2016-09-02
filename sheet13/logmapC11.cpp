#include <omp.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

// templated lyanpunov function returning type of c
template <typename T, typename S, typename W>
auto lyapunov(T c, S x, W l) -> decltype(c)
{
	auto sum = 0.;
	for (auto i = 0; i < l; ++i)
	{
		sum += log(std::abs(c*(1.-2.*x[i])));
	}
	return sum/l;
}

// Logistic map without explicit variable declarations
int main ()
{
	std::ofstream file1, file2;
	file1.open("chaos.txt", std::ios::trunc);
	file2.open("lyapunov.txt", std::ios::trunc);	

	auto len = 10000;

	// Auxiliary variables for making arrays
	auto aux = 0.;
	auto auxp = &aux;

	auto n = omp_get_num_procs();

	omp_set_num_threads(n);

	// auto arrays are not supported, but this is a hack to work around explicit variable declaration
	decltype(aux) c[len];

	// Coefficients
	for (auto i = 0; i < len; ++i)
	{
		c[i] = 1 + 3./9999.*i;
	}

	// Two more hacks
	auto A = new decltype(auxp)[len];
	auto L = new decltype(aux)[len];

	// Parallel section
	#pragma omp parallel shared(A,L)
	{
		// Parallel for loop
		#pragma omp for schedule(dynamic)
		for (auto i = 0; i < len; ++i)
		{
			// Thread-specific array
			auto values = new decltype(aux)[len];

			// Also, the undefined declaration is not supported for auto
			decltype(aux) x;
			
			
			for (auto j = 0; j < len; ++j)
			{
				// recursive for loop
				x = c[i]*x*(1.-x);
				values[j] = x;
			}
			// exporting values
			A[i] = values;
		}
		
		// calculating lyapunov exponents
		#pragma omp for schedule(dynamic)
		for (auto i = 0; i < len; ++i)
		{
			auto x = lyapunov(c[i],A[i],len);
			L[i] = x;
		}
		
	}

	// Print to file
	for (auto i = 0; i < len; ++i)
	{
		auto t = A[i];
		for (auto j = len-100; j < len; ++j)
			file1 << c[i] << "\t" << t[j] << std::endl;
		file2 << c[i] << "\t" << L[i] << std::endl;
	}

	file1.close();
	file2.close();	
}
