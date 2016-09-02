// Parallelized FFT
#include <fstream>
#include <iostream>
#include <complex>
#include "helfer.h"
#include <omp.h>
#include "fft.h"


// Heaviside functor class
class Functor
{
public:
	Functor() {}

	// the ()-operator, returns function value
	double operator () (const double x)
	{
		// returns sum of two return values: 1 or 0
		return ((x+45)>=0) - ((x-55)>= 0);
	}
};


int main ()
{
	// Sample limits
	double xmin = -100;
	double xmax = 100;
	
	// Files for printing results
	std::ofstream file_real, file_imag;
	file_real.open("Parallel-FFT-real.txt", std::ios::trunc);
	file_imag.open("Parallel-FFT-imag.txt", std::ios::trunc);
	
	// Number of parallel threads computing the FFTs
	int m = omp_get_num_procs();
	
	// Exponent
	int n = 24;

	// Number of samples
	int N = pow(2,n);	

	// Samples per thread	
	int M = N/m;

	// Matrix for storing FFT results
	Matrix<doz> T(M,m);
	
	// Stepsize
	double lambda = (xmax-xmin)/N;
	
	// Initialize threads
	omp_set_num_threads( m );

	// Begin parallel region
	#pragma omp parallel shared(xmin,xmax, T)
	{

		// Each thread has its own Heaviside function and vector for storing temporary results
		Functor f;
		VecDoz v(M,0.);
		int delta;

		// Thread number j < m
		int j = omp_get_thread_num();
		
		// Generating frequencies
		for (int J = 0; J < M; ++J)
		{
			delta = J*m+j;
			v[J] = f(xmin+lambda*delta);
		}

		// Fast Fourier Transform
		fft(v,1);

		// Multiplying element-wise with a phase factor
		for (int K = 0; K < v.size(); ++K)
		{
			double argument = 2.*M_PI*j*K/(M*m);	
			doz phase(cos(argument),sin(argument));
			v[K] *= phase;
		}

		// Storing results in matrix, one column per thread
		for (int i = 0; i < v.size(); ++i)
			T[i][j] = v[i];
		
		// Wait for all threads to finish, since T has to be fully formed for the next step
		# pragma omp barrier

		// Threaded for loop performing the last FFT for every row of the matrix
		# pragma omp for schedule(dynamic)
		for (int K = 0; K < M; ++K)
		{
			// Retrieving row
			VecDoz row(m,0.);
			for (int i = 0; i < m; ++i)
				row[i] = T[K][i];
			
			// Last Fast Fourier Transform and returning row to matrix
			fft(row,1);
			for (int i = 0; i < m; ++i)
				T[K][i]=row[i];
		}
	} // End parallel
		

	// Extracting real and imaginary values from the matrix
	VecDoub output_real(N,0.);
	VecDoub output_imag(N,0.);
	int index = 0;
	
	for (int k = 0; k < m; ++k)
	{
		for (int K = 0; K < M; ++K)
		{
			// Using the real and imag-functions provided by the <complex> library
			output_real[index] = real(T[K][k]);
			output_imag[index] = imag(T[K][k]);
			++index;
		}
	}

	// New parallel region. Two threads can write to file simultaneously.
	#pragma omp parallel num_threads(2)
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
	} // End parallel region

	file_real.close();
	file_imag.close();
}
