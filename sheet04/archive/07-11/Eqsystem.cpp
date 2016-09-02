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
	
	// Definition functor
	Functor func;
	
	// Two empty vectors for function values and jacobian
	double f[2];
	
	// initial values
	double x1[] = {1,1};
	double x2[] = {2,2};

	// function calls to mnewt() 
	int n1 = mnewt(func,x1,f);
	int n2 = mnewt(func,x2,f);
	
	// printing results
	cout << n1 << " " << x1[0] << " " << x1[1] << endl;
	cout << n2 << " " << x2[0] << " " << x2[1] << endl;
}
