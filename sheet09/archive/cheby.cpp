#include<cmath>
#include "helfer.h"


class Functor
{
public:
	double operator () (const double& x)
	{
		return x*sin(32.*x*x);
	}
};

class Cheby
{
	VecDoub coeffs;
	int m;
	
public:
	Cheby(const Functor&, const double&);
	Cheby(const VecDoub& V);
	
	double operator () (const double& x);
	
	Cheby quad();
};


Cheby::Cheby(const Functor& func, const double& acc=1e-10)
{
	
}

Cheby::Cheby(const VecDoub& V)
{
	
}

double Cheby::operator () (const double& x)
{
	double fval;
	for (int i = 0; i < m; ++i)
		...
	return fval;
}

Cheby Cheby::quad()
{
	VecDoub cint(n,0.);
	...
	return Cheby(cint);
}

