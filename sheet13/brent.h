#ifndef _BRENT_H_
#define _BRENT_H_

#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>

// Implementation of the van-Wijngaarden-Dekker-Brent method
// Numerical root finding of any functor that supports calls
// according to double func(const double x)

// **********
template <class T>
int brent(T& func, const double x1, const double x2, double& z, const double tol)
{
	// Maximum number of iterations
	static const int ITMAX = 100;
	static const double eps=std::numeric_limits<double>::epsilon();
	double a=x1, b=x2, c=x2, d, e;
	double fa=func(a), fb=func(b), fc, p, q, r, s, tol1, xm;

	if ((fa > 0. && fb > 0.) || (fa < 0. && fb < 0.))
	{
		std::cerr << "Root has not been bracketed!" << std::endl;
		return 0;
	}
	fc = fb;

	for (int i=0; i<ITMAX; ++i)
	{
		// Adjust the interval d
		if ((fb > 0. && fc > 0.) || (fb < 0. && fc < 0.))
		{
			c = a;
			fc = fa;
			e = d = b-a;
		}

		// Wrong direction: rename a, b, and c
		if (std::abs(fc) < std::abs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		// Convergency check
		tol1 = 2. * eps * std::abs(b) + 0.5*tol;
		xm = 0.5 * (c-b);

		if (std::abs(xm) <= tol1 || fb == 0.)
		{
			z = b;
			return i;
		}

		// Attempt quadratic interpolation
		if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
		{
			s = fb/fa;
			if (a == c)
			{
				p = 2. * xm * s;
				q = 1. - s;
			}
			else
			{
				q = fa/fc;
				r = fb/fc;
				p = s * ( 2. * xm * q * (q-r) - (b-a) * (r-1.) );
				q = (q-1.) * (r-1.) * (s-1.);
			}

			// Are we still inside the interval boundaries?
			if (p > 0.)
				q = -q;
			p = std::abs(p);

			double min1 = 3. * xm * q - std::abs(tol1*q);
			double min2 = std::abs(e*q);

			// Accept interpolation
			if (2.*p < (min1 < min2 ? min1 : min2))
			{
				e = d;
				d = p/q;
			}

			// Interpolation faile; try bisection
			else
			{
				d = xm;
				e = d;
			}
		}

		// Intervall shrinks too slowly; use bisection
		else
		{
			d = xm;
			e = d;
		}

		// Save last approximation to a
		a = b;
		fa = fb;

		// Test new value for the root
		if (std::abs(d) > tol1)
			b += d;
		else
			b += xm >= 0. ? (tol1 >= 0. ? tol1 : -tol1) : (tol1 >= 0. ? -tol1 : tol1);
		fb = func(b);
	}

	std::cerr << "Maximum number of iterations has been reached!" << std::endl;
	return 0;
}

#endif
