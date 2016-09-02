// moments.h
#ifndef _MOMENTS
#define _MOMENTS
#include "helfer.h"

class Moments
{
	int n;
	double av, va, sd, sk, cu;

	public:
// 	Constructor: calculate all moments of the data vector
	Moments(VecDoub& data);
	
	double ave();

	double sdev();

	double skew();

	double kurt();
};
#endif