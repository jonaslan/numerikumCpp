// moments.cpp

#include "moments.h"

Moments::Moments(VecDoub& data) : av(0.), va(0.), sd(0.), sk(0.), cu(0.)
{
	if ((n=data.size()) < 2)
		throw("Number of data must be at least 2!");

	for (int i=0; i<n; ++i)
		av += data[i];
	av /= n;

	double s, p, ep = 0.;
	for (int i=0; i<n; ++i)
	{
		ep += (s=data[i]-av);
		va += (p = s*s);
		sk += (p *= s);
		cu += (p *= s);
	}

	// Two-pass algorithm for the variance
	va = ( va - SQR(ep)/n )/(n-1);

	// Standard deviation
	sd = sqrt(va);

	if (va != 0.)
	{
		// Skewness
		sk /= (n * va * sd);

		// Kurtosis
		cu = cu/(n * va * va) - 3.;
	}
}
	
double Moments::ave()
{
	return av;
}

double Moments::sdev()
{
	return sd;
}

double Moments::skew()
{
	return sk;
}

double Moments::kurt()
{
	return cu;
}