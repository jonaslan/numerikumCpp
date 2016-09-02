#include <omp.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

template <typename T, typename S, typename W>
auto lyapunov(T c, S x, W l) -> decltype(c)
{
	auto sum = 0.;
	for (auto i = 0; i < l; ++i)
	{
		sum += log(std::abs(c*(1.-x[i])));
	}
	return sum/(l-1.);
}

int main ()
{
	std::ofstream file1, file2;
	file1.open("chaos.txt", std::ios::trunc);
	file2.open("lyapunov.txt", std::ios::trunc);	

	auto len = 10000;
	auto aux = 0.;
	auto auxp = &aux;

	auto n = omp_get_num_procs();

	omp_set_num_threads(n);

	decltype(aux) c[len];

	for (auto i = 0; i < len; ++i)
	{
		c[i] = 1 + 3./9999.*i;
	}

	auto A = new decltype(auxp)[len];
	auto L = new decltype(aux)[len];

	#pragma omp parallel shared(A)
	{
		#pragma omp for schedule(dynamic)
		for (auto i = 0; i < len; ++i)
		{
			auto values = new decltype(aux)[len];
			decltype(aux) x;			
			for (auto j = 0; j < len; ++j)
			{
				x = c[i]*x*(1.-x);
				values[j] = x;
			}
			A[i] = values;
		}
		
		#pragma omp for schedule(dynamic)
		for (auto i = 0; i < len; ++i)
		{
			auto x = lyapunov(c[i],A[i],len);
			L[i] = x;
		}
		
	}

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
