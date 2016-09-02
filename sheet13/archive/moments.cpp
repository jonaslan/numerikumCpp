// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 12 Problem 1 ***
// **************************

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "xorshift.h"
#include "helfer.h"
using namespace std;

// **********
class Moments
{
	int n;
	double av, va, sd, sk, cu;

	public:
// 	Constructor: calculate all moments of the data vector
	Moments(VecDoub& data) : av(0.), va(0.), sd(0.), sk(0.), cu(0.)
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

// 	Output functions
	double ave()
	{
		return av;
	}

	double sdev()
	{
		return sd;
	}

	double skew()
	{
		return sk;
	}

	double kurt()
	{
		return cu;
	}
};

// **********
int main()
{
// 	Initialize random seed value
	srand((unsigned)time(NULL));
	Xorshift rd(42);

	int i = 0, N = 1e6;

// 	Vector for the Gaussian random numbers
	VecDoub s(N);

// 	Polar generation method for Gaussian random numbers
	while (i<N/2)
	{
		double u1 = rd.doub(-1.,1.);
		double u2 = rd.doub(-1.,1.);

		double q = u1*u1 + u2*u2;

		double eps = numeric_limits<double>::epsilon();
		if (q <= eps || q >= 1.-eps)
			continue;

		double p = sqrt( -2. * log(q)/q );

		s[2*i]   = p * u1;
		s[2*i+1] = p * u2;
		++i;
	}

// 	Initialize object for the calculation of the moments
	Moments m(s);
	cout << "Mean value: " << m.ave() << endl;
	cout << "Std. dev.:  " << m.sdev() << endl;
	cout << "Skewness:   " << m.skew() << endl;
	cout << "Kurtosis:   " << m.kurt() << endl;

// 	Initialize bins for the distribution function
	int M = 256;
	VecInt bin(M,0);
	const double a = 3.;
	double dx = 2.*a/M;

// 	Sort particles by calculating the correct bin index
	for (int i=0; i<N; ++i)
		if (abs(s[i]) < a)
			++bin[ floor(M/2+s[i]/dx) ];

// 	Output the resulting distribution to a file
	ofstream datei("moments.txt", ios::trunc);

	double x = -a+0.5*dx;
	for (int i=0; i<M; ++i, x+=dx)
		datei << x << "\t" << bin[M/2]*exp(-x*x/2.) << "\t" << bin[i] << endl;

	datei.close();
}
