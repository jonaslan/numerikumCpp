#include <iostream>
#include <cstdlib>
#include <cmath>

// golden section

using namespace std;


// function f(x) = (x-1)²+x
double f(const double& x)
{
	return (x-1)*(x-1) + x;
}


int minimize(double& a, double& b, double& c)
{
	static const double G = (-1.+sqrt(5.))/2.;
	double x,bc,ab;
	double diff = 1;
	int count = 0;

	// loop with specified break conditions
	while (diff > 10e-8 && count <= 1000)
	{
		// distances between 
		bc = abs(b-c);
		ab = abs(a-b);
		
		// if b is closer to c than a:
		if (bc > ab)
		{
			// place x between c and b
			x = b + (1-G)*(c-b);
		} else {
			// else between a and b
			x = a + G*(b-a);
		}
		
		// if b is smaller than x:
		if (b < x)
		{
			// and the function value of x is larger than that of b, i.e
			// b is probably closer to the minimum of the function:
			if (f(x) > f(b))
			{
				// move the point x closer to the others
				c = x;
				diff = abs(b-x);
			} else {
			// else, if f(x) < f(b) -> x probably closer to minimum ->
			// move a closer to the other points
				a = b;
				diff = abs(b-x);
				b = x;
			}
		} else {
		// equivalent but opposite if x > b
			if (f(x) > f(b))
			{
				a = x;
				diff = abs(b-x);
			} else {
				c = b;
				diff = abs(b-x);	
				b = x;
			}
		}
		++count;
	}
	return count;
}

int main ()
{
	//starting values
	double a = -4.;
	double b = 3.;
	double c = 9.;
	
	// function call and output
	cout << "Count: " << minimize(a,b,c) << endl;
	cout << "a = " << a << " b = " << b << " c = " << c << endl;
	cout << "Function value at minimum: " << f(a) << endl;
	
}