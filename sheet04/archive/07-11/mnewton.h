#ifndef _MNEWTON_H_
#define _MNEWTON_H_

#include <cmath>

// Implementation of the multidimensional Newton-Raphson method
// Numerical root finding of any functor that supports calls
// according to func(double x[], double f[], double p[])

// **********
template <class T>
int mnewt(T& func, double x[], double f[])
{
	static const int n = 2;
	static const int ntrial = 1000;
	static double tolx = 1e-15;
	static double tolf = 1e-15;

	double p[n];
	for (int i=0; i<ntrial; ++i)
	{
		func(x,f,p);

		double errf = 0.;
		for (int j=0; j<n; ++j)
			errf += std::abs(f[j]);

		if (errf <= tolf)
			return i;

		double errx = 0.;
		for (int j=0; j<n; ++j)
		{
			errx += std::abs(p[j]);
			x[j] += p[j];
		}

		if (errx <= tolx)
			return i;
	}

	std::cerr << "Maximum number of iterations has been reached!" << std::endl;
	return 0;
}

#endif
