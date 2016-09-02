		
// FFT of rectangular function

#include <iostream>
#include <fstream>
#include "fft.h"

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
	// boundaries for x
	const double xmin = -100;
	const double xmax = 100;

	Functor f;
	std::ofstream file1, file2;
	file1.open("heaviside-real.txt", std::ios::trunc);
	file2.open("heaviside-imaginary.txt", std::ios::trunc);
	
	// the number of steps in the interval is 2^n
	int n = 10;
	int N = pow(2,n);
	
	// stepsize
	double lambda = (xmax-xmin)/N;
	
	VecDoub data(2*N,0.);	
	
	// fill data VecDoub with function values
	for (int i = 0; i < 2*N; i+=2)
	{
		data[i] = f(xmin+lambda*i);
	}
	
	fft(data,1);

	 //Print to file	 
	
	 double freq = 0;
	 for (int i = 0; i < 2*N; i +=2)
	 {
		// frequency increases to nyquist frequency, i.e. 1/(2*lambda)
		if (i >= 2 && i <= N)
		{
			freq = freq+1./(N*lambda);
		}

		// frequency changes sign at nyquist frequency
		if(i == N)
			freq = -freq;

		// and then starts decreasing until 0
		else if (i > N)
		{
			freq = freq + 1./(N*lambda);
		}
		
		file1 << freq << "\t" << data[i] << std::endl;
		file2 << freq << "\t" << data[i+1] << std::endl;
	}
	
	
	
	fft(data,-1);	

	// checking back transformation quality (the data values are corrected for discretization	
	//for (int i = 0; i < 2*N; i+=2)
	//	std::cout << data[i]/N << "\t" << f(xmin+lambda*i) << std::endl;
	file1.close();
	file2.close();

	
}
