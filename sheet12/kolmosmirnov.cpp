#include <iostream>
#include <cmath>
#include "helfer.h"


// Exponential distribution function with expected value 100
class Functor 
{
public:
	Functor() {};
	
	double operator () (const double x)
	{
		return 1. - exp(-x/100.);
	}
};


// the P function
class P
{
public:
	P() {};

	double operator () (const double z)
	{
		double sum = 0;
		// the sum converges extremely fast
		for (int j = 1; j < 100; ++j)
		{
			sum += exp((-(2.*j-1.)*(2.*j-1.)*M_PI*M_PI)/(8.*z*z));
		}
		return std::sqrt(2.*M_PI)/z*sum;
	}
};

int main ()
{
	Functor F;
	P p;

	// sample values
	const double X[] = {10., 72., 81., 94., 112., 116., 124., 140., 145., 155.};

	// number of samples
	const double N = 10.;

	// the two auxiliary quantities
	double dl = 0.;
	double dr = 0.;


	// Finding max values for the deviation from the analytic distribution and assigning these values to dr & dl
	for (int i = 0; i < 10; ++i)
	{		
		double d1 = (i+1)/N - F(X[i]);
		if (d1 > dl)
			dl = d1;
		double d2 = F(X[i]) - (i)/N;
		if (d2 > dr)
			dr = d2;
	}
	
	// the Dks-factor is the maximum value of dr & dl
	double Dks = MAX(dr,dl);

	// full z-value
	double z = ( std::sqrt(N) + 0.12 + 0.11/std::sqrt(N) )*Dks;

	// Printing probability
	std::cout << "The probability of observing a deviation larger than " << Dks << " is: "<< 1-p(z) << std::endl;
}
