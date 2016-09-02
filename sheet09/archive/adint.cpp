#include <cmath>
#include <iostream>
#include <fstream>
#include "adapt.h"

class Functor
{
	double lambda, s, t,eps;
	int count;

double f(const double x)
{
    if ( x > sqrt(3.)/(2.*lambda)+eps )
    {
        return (pow((1+x*x),-s/2.)/sqrt(-1./(4.*lambda*lambda)+x*x/3.))
		*sin(t*sqrt(-1./(4.*lambda*lambda)+x*x/3))*exp(-t/(2*lambda));
    }
    else if ( x < sqrt(3.)/(2.*lambda)-eps )
    {
        return (pow((1+x*x),-s/2.)/sqrt(1./(4.*lambda*lambda)-x*x/3.))
		*sinh(t*sqrt(1./(4.*lambda*lambda)-x*x/3))*exp(-t/(2*lambda));
    }
    else
    {
        return x;
    }


}

public:

	Functor() : lambda(3.), s(5./3.), t(0), eps(1e-10), count(0) {}

	double operator () (double k)
	{
        if(k == 0)
            k = 1e-20;
        return f(k)+f(1./k)*(1./k)*(1./k);
	}

	void setTime(const double tin)
	{
		t = tin;
	}


};


int main ()
{
	Functor f;
	Adapt A(1e-8);
	ofstream file;

	file.open("Adapt.txt");

	for (int i = 1; i <= 100; ++i)
	{
		f.setTime((double)i);
		file << i << "\t" << A.integrate(f,0.,1.0) << std::endl;
	}
}

