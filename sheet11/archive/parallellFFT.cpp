// Parallelized FFT
#include <fstream>
#include <iostream>
#include <complex>
// #include <vector>
#include "helfer.h"
#include <omp.h>
#include <ctime>
#include "fft.h"

// OMP Commands:
// ifdef _OPENMP
// omp_set_num_threads( omp_get_num_procs() );
// pragma omp parallel shared(var1,var2...)
// 
// pragma omp for schedule(schedule)
// 
// pragma omp critical
// 
// 

class Functor
{
public:
	Functor() {}

	// the ()-operator, returns function value
	double operator () (const double x)
	{
		return ((x+45)>=0) - ((x-55)>= 0);
	}
};


int main ()
{
	double xmin = -100;
	double xmax = 100;
	
	std::ofstream file_real, file_imag;
	file_real.open("Parallel-FFT-real.txt", std::ios::trunc);
	file_imag.open("Parallel-FFT-imag.txt", std::ios::trunc);
	
	int m = omp_get_num_procs();
	int n = 10;
	int N = pow(2,n);	
	int M = N/m;
	Matrix<doz> T(M,m);
	
	double lambda = (xmax-xmin)/N;
	
	omp_set_num_threads( m );
	#pragma omp parallel shared(xmin,xmax, T)
	{
		Functor f;
		
		VecDoz v(M,0.);
		int delta;
		int j = omp_get_thread_num();
		
		for (int J = 0; J < M; ++J)
		{
			delta = J*m+j;
			v[J] = f(xmin+lambda*delta);
		}
		fft(v,1);
		for (int K = 0; K < v.size(); ++K)
		{
			double argument = 2.*M_PI*j*K/(M*m);	
			doz phase(cos(argument),sin(argument));
			v[K] *= phase;
		}

		for (int i = 0; i < v.size(); ++i)
			T[i][j] = v[i];
		
		# pragma omp barrier
		# pragma omp for schedule(dynamic)
		for (int K = 0; K < M; ++K)
		{
			VecDoz row(m,0.);
			for (int i = 0; i < m; ++i)
				row[i] = T[K][i];
			
			fft(row,1);
			for (int i = 0; i < m; ++i)
				T[K][i]=row[i];
		}
	} // End parallel
	
	
	VecDoub output_real(N,0.);
	VecDoub output_imag(N,0.);
	int index = 0;
	
	clock_t start,end;
	
	start = clock();
	for (int k = 0; k < m; ++k)
	{
		for (int K = 0; K < M; ++K)
		{
			output_real[index] = real(T[K][k]);
			output_imag[index] = imag(T[K][k]);
			++index;
		}
	}
	end = clock();
	
	std::cout << (end-start)/double(CLOCKS_PER_SEC) << std::endl;
	
	#pragma omp parallel
	{
		double freq = 0;
		#pragma omp for schedule(dynamic)
		for (int i = 0; i < N; ++i)
		{
// 			frequency increases to nyquist frequency, i.e. 1/(2*lambda)
			if (i >= 1 && i <= N/2)
			{
				freq = (1./(N*lambda))*i;
			}

// 				and then starts decreasing until 0
			else if (i > N/2)
			{
				freq = -1/(2*lambda) + (i-(N/2))*(1./(N*lambda));
			}
			#pragma omp critical
			file_real << freq << "\t" << output_real[i] << std::endl;
			#pragma omp critical
			file_imag << freq << "\t" << output_imag[i] << std::endl;
		}
	}
	file_real.close();
	file_imag.close();
}
