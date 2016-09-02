// Bracketing and bisection

// N = 100, c = 1,1


#include <iostream>
#include <cmath> //exp()
#include<cstdio> // printf()

using namespace std;

// global variables
const static int N = 100;
const static double redFactor = 1.1;

// Functor class declaration and definition
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

// bracketing function, takes functor and boundary values by reference
int bracket (Functor& func, double& a, double& b)
{
	// function values for boundaries
	double f = func(a);
	double g = func(b);
	
	// a maximum of 100 iterations
	for (int i = 0; i < N; ++i)
	{
		
		// are the function values on opposite sides of 0 -> bracketed
		if (f*g < 0)
		{
			return i;
		}
		// if function value of a is smaller than that of b (possibly closer to zero)
		// -> move a "downhill"; to the right if a is bigger than b, else to the left,
		// and take new function value
		if (abs(f) < abs(g))
		{
			f = func(a += redFactor*(a-b));
		}
		// else do same procedure for b -> new function value g(b).
		else
		{
			g = func(b += redFactor*(b-a));
		}
	}
	return 0;
}

// bisection method root finding function
int root(Functor& myfunc, double& a, double&b, double& ans)
{
		double mid = (a+b)/2.;
		double diff = mid-a;
		int count = 0;
	
		// while loop restricted to accuracy and runs
		while (diff > 1e-14 && count < 1000)
		{
			// bracket function call, checking if a and midpoint bracket zero
			if (bracket(myfunc,a,mid))
			{
				// set b to midpoint, take new midpoint and calculate difference
				b = mid;
				mid = (a+b)/2.;
				diff = mid-a;
			}
			// if a and midpoint 
			else
			{
				a = mid;
				mid = (a+b)/2.;
				diff = mid-a;
			}
			++count;
		}
		// set answer variable to midpoint in bisection
		ans = mid;
	
		// return number of iterations
		return count;
}

int main () {
	Functor myfunc(4.);
	double a = 1.0;
	double b = 1.5;
	double zero;
	
	int iter = bracket(myfunc, a, b);
	
	if (iter == 0)
	{
		cout << "Failed to bracket\n";
	}
	else
	{
		cout << "Bracketed in " << iter << " iterations by using: ";
		cout << "a: " << a << " b: " << b << endl;
	
		// bisection
	
	iter = root(myfunc, a, b, zero);
	double fval = myfunc(zero);

	printf("Root: %f, Function value at root: %f, reached in: %d iterations", zero, fval, iter);
	}
}
