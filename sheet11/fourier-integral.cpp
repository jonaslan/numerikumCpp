#include <iostream>
#include <fstream>
#include <cmath>
#include "dftintegrate.h"
#include "helfer.h"

// Definition of the class for the function we want to integrate with the Fourier kernel
class Functor
{
private:
	double a;

public:
	Functor() : a(2) {}

	double operator () (double x)
	{
		return exp(-a*x*x); 
	}
};

int main () {
	// Lower and upper limits of integration
	const double a = -10;
	const double b = 10;
	Functor f; // Instantiate the functor
	
	// Open file
	std::ofstream file1;
	file1.open("fourier-integral-approximation.dat", std::ios::trunc);
	// Instantiation of Fourier integral. Note that the constructor performs the Fourier transform calculations
	FFTint<Functor> h(f,a,b);

	double k_min = 0, k_max = 10; // Range of frequencies we are interested in
	const int k_steps = 100; // Steps in k over [0,10] for interpolation and plotting
	const double dk = (k_max - k_min)/k_steps;

	for (int i = 0; i <= k_steps; ++i)
	{
		double k = k_min + i*dk;
		double gr, gi; //Initialise real and imaginary parts of g(k)
		h.interp(k, gr, gi); // Interpolation
		double analytical_result = sqrt(M_PI/2)*exp(-k*k/8); // Calculation of analytical result (Gaussian function)
		file1 << k << '\t' << gr << '\t' << gi << '\t' << analytical_result << endl; // Write results to file
	}

	file1.close();
}
