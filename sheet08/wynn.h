#ifndef _WYNN_H_
#define _WYNN_H_

#include "helfer.h"

// **********
class Wynn
{
	VecDoub e;
	int n, ncv;
	double eps, small, big, lastval, lasteps;

	public:
	bool conv;

// 	Constructor: call with maximum number of terms, nmax,
// 	and with desired accuracy, epss
	Wynn(int nmax, double epss) : e(nmax), n(0), ncv(0), conv(false), eps(epss), lastval(0.)
	{
		small = numeric_limits<double>::min() * 10.;
		big = numeric_limits<double>::max();
	}

// 	Parameter: subsequent partial sums
// 	Return value: estimate for the true value of the sum
// 	If convergence if found, conv is set to true
	double next(double sum)
	{
		double diff, tmp1, tmp2, val;
		e[n] = sum;
		tmp2 = 0.;

		for (int j=n; j>0; --j)
		{
			tmp1 = tmp2;
			tmp2 = e[j-1];
			diff = e[j] - tmp2;

			if (abs(diff) <= small)
				e[j-1] = big;
			else
				e[j-1] = tmp1 + 1./diff;
		}

		++n;
		val = (n&1) ? e[0] : e[1];

		if (abs(val) > 0.01 * big)
			val = lastval;

		lasteps = abs(val-lastval);

		if (lasteps > eps)
			ncv = 0;
		else
			++ncv;

		if (ncv >= 3)
			conv = true;

		return (lastval = val);
	}
};

#endif
