#include "adapt.h"
#include "brent.h"
//#include <iostream>


int main ()
{
	const double min = 0.;
	const double max = 2.;
	Adapt A(1e-15);
	
	
	// the A.integrate function takes its first input by reference. Rvalues cannot usually be
	// taken by reference since they can't be changed.
	A.integrate([](double x) {double z = x*x; return z*z*std::log(x+std::sqrt(z+1));}, max, min);
	
	
	// C++11:
	
	double arr[] = {3.,3.5,4.,4.5,5.,5.5,6.};
	double zero = 0;
	
	for (auto&& a : arr)
	{
		// Here the anonymous function is named, and can after that be referenced
		auto func = [=](double x) {return exp(x)-a*x;};
		brent(func,min,max,zero,1e-15);
	}
}