#include <iostream>
#include <cmath>
#include "mnewton.h"

using namespace std;


// Definition of functor class
class Functor
{
	public:	
	void operator ()(double x[], double f[], double p[])
	{
		// F(X)
	
		f[0] = 1 - cos(x[0]) - cos(x[1]);
		f[1] = 2*x[0] - x[1];
	
		// -J⁻1*F(X) == P
	
		p[0] = -(cos(x[0])+ cos(x[1]) - 1 + x[1]*sin(x[1]) - 2*x[0]*sin(x[1]))/(-sin(x[0])-2*sin(x[1]));
		p[1] = -(2*cos(x[0]) + 2*cos(x[1]) - 2 + 2*x[0]*sin(x[0]) - x[1]*sin(x[0]))/(-sin(x[0])-2*sin(x[1]));
	}
};


int main () {
	
	// Definition of call parameters
	Functor func;
	double a = 0.2379;
	double b = 1.2;
	double zero;
	double acc = 1e-14;
	
	double f[2];
	double p[2];
	double x[] = {1,1};
	
	func(x,f,p);

	int n = mnewt(func,x,f);
	
	cout << n << " " << x[0] << " " << x[1] << endl;
}
