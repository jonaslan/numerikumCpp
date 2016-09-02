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
	Quadrature() : n(0)
	{
	}

// 	Abstract method to be redefined in any class derived from Quadrature
	virtual double next() = 0;
};

// **********
template <class T>
class Trapzd : public Quadrature
{
	
	double out;
	double h,a,b;
	T& f;
public:

	Trapzd(T& func, const double ain, const double bin, int N=1) : out(0), a(ain), b(bin), h(b-a), f(func) {}

// Trap.next() calculates new intermediate points and evaluates the function at them, yielding a higher-order approximation to the integral	
	double next()
	{
		++n;
// Special case n = 1
		if (n == 1)
		{
			out = 0.5*(b-a)*( f(a) + f(b) );
			return out;
		}
	// implementation 2^(n-2) factor -> necessary to calculate step size
		int i = 1 << (n-2);
		
		h = (b-a)/i; // step size -> decreases as order increases
		
		double x = a + 0.5*h; // calculation of new intermediate points
		
		double sum = 0;
// Calculation of sum of function values at the new points, which is afterwards used (together with the previous estimate) to calculate the current estimate
		for (int j = 0; j < i; ++j, x+=h)
			sum += f(x);
		out = 0.5*( out + h*sum );
		return out;
	}
};

// The following function qtrap is integrates f with an increasing number of steps using Trapezoidal Rule with the stopping criteria:
// relative difference between successive approximations < eps or iteration > certain number of iterations

template <typename T>
void qtrap(T& f, const double a, const double b, int& steps, double& out, const double eps=1e-10)
{
	Trapzd<T> Trap(f,a,b);
	double integral,old,e;
	int count = 0;
	old = ~0;
	e = old;
// The loop contains successive calls to Trap.next providing higher order approximations	
	while (abs(e) > eps && count < 50)
	{
		std::cout.precision(20);
		std::cout << count << " " << integral << std::endl;
		integral = Trap.next();
		
		
		// relative error to check for convergence
		e = (old - integral)/old;
		
		old = integral;
		
		++count;
	}
	
	steps = count;
	out = integral;
}

// qsimp implements Simpson's rule (based on the trapezoidal rule)
template <typename T>
void qsimp(T& f, const double a, const double b, int& steps, double& out, const double eps=1e-10)
{
	Trapzd<T> Trap(f,a,b);
	
	double trapcurrn, trapcurr2n, trapoldn, trapold2n, oldsi, newsi, e;
	e = 10;
	trapold2n = 0;
	trapoldn = 0;
	int count = 0;
	
	while (abs(e) > eps && count < 50)
	{

		trapcurrn = Trap.next();
		trapcurr2n = Trap.next();
		
		newsi = ((4./3.)*trapcurr2n-(1./3.)*trapcurrn);
		oldsi = ((4./3.)*trapold2n-(1./3.)*trapoldn);
		
		e = (newsi-oldsi)/oldsi;
		
		std::cout.precision(20);
		std::cout << count << " " << newsi << std::endl;
		
		trapold2n = trapcurr2n;
		trapoldn = trapcurrn;
		++count;
	}
	
	out = newsi;
	steps = count;
}
	

#endif
