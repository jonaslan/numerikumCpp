// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 13 Problem 3 ***
// **************************

#ifndef _VDP_H_
#define _VDP_H_

#include "helfer.h"

// **********
class VanDerPol
{
	const double eps;

	public:
	VecDoub ystart;

	VanDerPol(const double epsIn = 1e-2);
	void operator() (const double x, const VecDoub& y, VecDoub& dydx);
};

#endif
