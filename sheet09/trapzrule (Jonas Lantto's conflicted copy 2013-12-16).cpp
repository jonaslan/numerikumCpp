#include <iostream>
#include <cmath>
#include "quadrature.h"

// The functors, one for each function
class Functor1
{
public:
	
	double operator() (const double x)
	{
		return x*(1+x*x);
	}
};

class Functor2
{
public:
	
	double operator() (const double x)
	{
		return x*x*x*x*log(x + std::sqrt(x*x + 1));
	}
};

int main ()
{
// declarations, range, error, output variables
	Functor1 f1;
	Functor2 f2;
	double a,b,eps;
	double out1,out2;
	int steps1,steps2;
	eps = 1e-8;
	a = 0;
	b = 2;
	
	// quadrature trapezoid rule function calls
	qtrap<Functor1>(f1,a,b,steps1,out1,eps);
	qtrap<Functor2>(f2,a,b,steps2,out2,eps);
	
	std::cout << out1 << " steps: " << steps1 << std::endl << out2 << " steps: " << steps2 << std::endl;
	
	// Simpson's method functions calls
	qsimp<Functor1>(f1,a,b,steps1,out1,eps);
	qsimp<Functor2>(f2,a,b,steps2,out2,eps);
	
	std::cout << out1 << " steps: " << steps1 << std::endl << out2 << " steps: " << steps2 << std::endl;
	
}