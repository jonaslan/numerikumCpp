#include <iostream>
#include "brent.h"

using namespace std;


// Definition of functor class
class Functor
{
	double a;

	public:
	
	Functor(double ain)
	{
		a = ain;
	}
	double operator () (double x)
	{
		return exp(x)-a*x;
	}
};

int main () {
	
	// Definition of call parameters
	Functor func(4.);
	const double a = 0.2379;
	const double b = 1.2;
	double zero;
	double acc = 1e-14;
	
	// Function call
	int n = brent(func,a,b,zero,acc);
	cout << "Zero: " << zero << " in " << n << " iterations.\n";
	
	// Brent's method: 6 iterations
	// Bisection method: 135 for the same accuracy level
	
}
