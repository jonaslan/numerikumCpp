#ifndef _PADE_H_
#define _PADE_H_

#include "helfer.h"
#include "ludcmp.h"
using namespace std;

// **********
class Pade
{
	int N;
	VecDoub a, b;

	public:
// 	Constructor: calculate Pade coefficients based on Taylor coefficients koef
	Pade(const VecDoub& koef) : N((koef.size()-1)/2)
	{
		a.resize(N+1);
		b.resize(N+1);

		VecDoub x(N), y(N);
		MatDoub c(N,N);

		for (int k=0; k<N; ++k)
		{
			y[k] = -koef[N+k+1];
			for (int m=0; m<N; ++m)
				c[k][m] = koef[N-m+k];
		}

		// Matrix inversion
		LUdcmp lu(c);
		lu.solve(y,x);

		a[0] = koef[0];
		b[0] = 1.;

		// Determine coefficient arrays a and b
		for (int k=0; k<N; ++k)
		{
			b[k+1] = x[k];
			double sum = koef[k+1];
			for (int j=0; j<=k; ++j)
				sum += x[j] * koef[k-j];
			a[k+1] = sum;
		}
	}

	// ()-operator overload
	double operator() (const double& x)
	{
		// same procedure as in the Taylor approximation with multiplication
		// by x and summation by the next term in the koef-VecDoub
		double suma = a[N];
		double sumb = b[N];
		for (int i = N-1; i >= 0; --i)
		{
			suma *= x;
			suma += a[i];
			sumb *= x;
			sumb += b[i];
		}
		// return quotient of the sums (R(x))
		return suma/sumb;
	}

};

#endif
