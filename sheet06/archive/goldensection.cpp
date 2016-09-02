#include <iostream>
#include <cstdlib>
#include <cmath>

// golden section

using namespace std;

double f(double& x)
{
	return (x-1)*(x-1) + x;
}

static double G = (-1.+sqrt(5.))/2.;

int minimize(double& a, double& b, double& c)
{

	double x,bc,ab;
	double diff = 1;
	int count = 0;
	
	while (diff > 10e-8 && count <= 1000)
	{
		bc = abs(b-c);
		ab = abs(a-b);
		if (bc > ab)
		{
			x = c - G*(c-b);
		} else {
			x = a + G*(b-a);
		}
		if (b < x)
		{
			if (f(x) > f(b))
			{
				c = x;
			} else {
				a = b;
				diff = abs(b-x);
				b = x;
			}
		} else {
			if (f(x) > f(b))
			{
				a = x;
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
	
	cout << "Count: " << minimize(a,b,c) << endl;
	cout << "a = " << a << " b = " << b << " c = " << c << endl;
	
}