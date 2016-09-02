#include <iostream>
#include "helfer.h"
#include<cmath>
#include <fstream>
#include "pade.h"

using namespace std;

class Taylor {

private:

	VecDoub coeffs;
	
public:
	
	Taylor(const VecDoub& t) : coeffs(t) {}
	
	double operator () (const double& x)
	{
		// iterate from the back of the coefficients and multiply by x and add
		// the next coefficient in every iteration. 
		double returnval = coeffs[4];
		for (int i = 3; i >= 0; --i)
		{
			returnval *= x;
			returnval += coeffs[i];
		}
		return returnval;
	}
	// The Taylor approximation diverges from the analytic solution for x-values 
	// above ~2. For x < 1 the absolute error is under 0.001. 
};


// f(x) - the analytic expression
double function(const double& x)
{
	return pow(pow(1.+x,4./3.)+7.,1./3.);
}

int main(){

	std::ofstream file;
	file.open("pade-plot.txt",std::ios::trunc);
	
	// a0 = 2
	double a[] = {2., 1./9., 1./81., -49./8748., 175./78732.};
	
	VecDoub c(5,a);
	Taylor T(c);
	Pade P(c);
	
	double x = 0.;
	
	for (int i = 0; i <= 500; ++i)
	{
		double f = function(x);
		double t = T(x);
		double p = P(x);
		
		file << x << "\t" << f << "\t" << t << "\t" << p << std::endl;
		x += 0.02;
	}
	file.close();	
	
	// The plot shows that the Pade approximant retains its validity over a 
	// broader range of values than the Taylor approximation
	
	
	
	
}