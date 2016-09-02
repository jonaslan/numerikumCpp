#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_
#include <cmath>
#include <iostream>

using namespace std;

// **********
// Abstract base class Quadrature
class Quadrature
{
	protected:
	int n;

	public:
	Quadrature() : n(0)	{}

// 	Abstract method to be redefined in any class derived from Quadrature
	virtual double next() = 0;
};

// **********
template <class T>
class Trapzd : public Quadrature
{
	
	// private variables and template pointer to functor
	double out;
	double h,a,b;
	T& f;
public:

	// constructor with initialization list
	Trapzd(T& func, const double ain, const double bin, int N=1) : out(0), a(ain), b(bin), h(b-a), f(func) {}

	// next function
	double next()
	{
		// begin by incrementing n
		++n;
		// if n is one return average of the start and end values of the function
		if (n == 1)
		{
			out = 0.5*(b-a)*( f(a) + f(b) );
			return out;
		}
		
		// multiply number of steps by 2 for every level level of accuracy
		int i = 1 << (n-2);
		
		// calculate step size
		h = (b-a)/i;
		
		// first x-value in range
		double x = a + 0.5*h;
		
		// sum up all function values for all steps in interval
		double sum = 0;
		for (int j = 0; j < i; ++j, x+=h)
			sum += f(x);
		// add sum to previous sum and return the average
		out = 0.5*( out + h*sum );
		return out;
	}
};

// template quadrature trapezoid rule implementation function
template <typename T>

// function takes template function, range values, iteration variable, output variable and precision value
void qtrap(T& f, const double a, const double b, int& steps, double& out, const double eps=1e-10)
{
	// instantiation of trapezoid solver class
	Trapzd<T> Trap(f,a,b);
	
	// variables
	double integral,old,e;
	int count = 0;
	old = 1e10;
	e = old;
	
	// while loop for calculating integral 
	while (abs(e) > eps && count < 50)
	{
		std::cout.precision(20);
		std::cout << count << " " << integral << std::endl;
		
		// new and improved integral
		integral = Trap.next();
		
		// relative error
		e = (old - integral)/old;
		
		old = integral;
		
		++count;
	}
	
	// set output variables
	steps = count;
	out = integral;
}

// template Simpson function
template <typename T>
void qsimp(T& f, const double a, const double b, int& steps, double& out, const double eps=1e-10)
{
	Trapzd<T> Trap(f,a,b);
	
	// the integral value holders for Simpson's rule
	double trapcurrn, trapcurr2n, trapoldn, trapold2n, oldsi, newsi, e;
	e = 10;
	trapold2n = 0;
	trapoldn = 0;
	int count = 0;
	
	while (abs(e) > eps && count < 50)
	{
		
		// new (current) trapezoid values
		trapcurrn = Trap.next();
		trapcurr2n = Trap.next();
		
		// Simpson approximation
		newsi = ((4./3.)*trapcurr2n-(1./3.)*trapcurrn);
		oldsi = ((4./3.)*trapold2n-(1./3.)*trapoldn);
		
		// relative error
		e = (newsi-oldsi)/oldsi;
		
		std::cout.precision(20);
		std::cout << count << " " << newsi << std::endl;
		
		// set old integral values
		trapold2n = trapcurr2n;
		trapoldn = trapcurrn;
		++count;
	}
	
	// set output variables
	out = newsi;
	steps = count;
}
	

#endif
