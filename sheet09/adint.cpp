#include <cmath>
#include <iostream>
#include <fstream>
#include "adapt.h"

// functor class
class Functor
{
	double lambda, s, t,eps;
	int count;

	
	// internal function 
	double f(const double x)
	{
		// if input value is larger than sqrt(3)/(2*lambda) return circular sinus version of function
		if ( x > sqrt(3.)/(2.*lambda)+eps )
		{
			return (pow((1+x*x),-s/2.)/sqrt(-1./(4.*lambda*lambda)+x*x/3.))
			*sin(t*sqrt(-1./(4.*lambda*lambda)+x*x/3))*exp(-t/(2*lambda));
		}
		// else, if it's smaller, return the hyperbolic sinus version.
		else if ( x < sqrt(3.)/(2.*lambda)-eps )
		{
			return (pow((1+x*x),-s/2.)/sqrt(1./(4.*lambda*lambda)-x*x/3.))
			*sinh(t*sqrt(1./(4.*lambda*lambda)-x*x/3))*exp(-t/(2*lambda));
		}
		// else, if the value is close to the singularity value, return the limit: x
		else
		{
			return x;
		}
	}

public:

	// constructor, sets constant and initial values
	Functor() : lambda(3.), s(5./3.), t(0), eps(1e-10), count(0) {}

	double operator () (double k)
	{
		// if k is zero, there will be a divide by zero error.
        if(k == 0)
            k = 1e-20;
			
		// return improper integral reformulation
        return f(k)+f(1./k)*(1./k)*(1./k);
	}

	// set time function
	void setTime(const double tin)
	{
		t = tin;
	}


};


int main ()
{
	// instantiations
	Functor f;
	// 1e-8 is moreover the practical limit for calculating the integral in a reasonable time
	Adapt A(1e-8);
	ofstream file;

	file.open("Adapt.txt");

	for (int i = 1; i <= 100; ++i)
	{
		// set new time and calculate integral
		f.setTime((double)i);
		file << i << "\t" << A.integrate(f,0.,1.0) << std::endl;
	}
}

