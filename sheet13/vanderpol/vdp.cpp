// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 13 Problem 3 ***
// **************************

#include "vdp.h"

// **********
VanDerPol::VanDerPol(const double epsIn) : eps(epsIn), ystart(2)
{
	ystart[0] = 2.;
	ystart[1] = 0.;
}

// **********
void VanDerPol::operator() (const double x, const VecDoub& y, VecDoub& dydx)
{
// 	van der Pol's differential equation
	dydx[0] = y[1];
	dydx[1] = ( 1.-SQR(y[0]) ) * y[1] - y[0];
	dydx[1] /= eps;
}
