#include <iostream>
#include <cmath>
#include "quadrature.h"
// Definition of functors
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
	//Instantiation of functors
	Functor1 f1;
	Functor2 f2;
	double a,b,eps;
	//Declaration of variables that will store final result of the integration
	double out1,out2;
	int steps1,steps2;
	eps = 1e-8; // Accuracy demanded
	// Limits of integration
	a = 0;
	b = 2;

	// Calculation of integral approximations by Trapezoidal Rule
	qtrap<Functor1>(f1,a,b,steps1,out1,eps);
	qtrap<Functor2>(f2,a,b,steps2,out2,eps);
	cout << "Trapezoidal rule:" << endl;	
	cout << "f1(x)\tIntegral result: " << out1 << " steps: " << steps1 << std::endl; 
	cout << "f2(x)\tIntegral result: " << out2 << " steps: " << steps2 << std::endl;

	// Calculation of integral approximations by Simpson's Rule
	qsimp<Functor1>(f1,a,b,steps1,out1,eps);
	qsimp<Functor2>(f2,a,b,steps2,out2,eps);


	cout << "Simpson's rule:" << endl;
	cout << "f1(x)\tIntegral result: " << out1 << " steps: " << steps1 << std::endl; 
	cout << "f2(x)\tIntegral result: " << out2 << " steps: " << steps2 << std::endl;

}
