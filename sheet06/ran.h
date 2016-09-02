// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 06 Problem 4 ***
// **************************

#ifndef _RAN_H_
#define _RAN_H_

#include <ctime>
#include <cstdlib>
#include "helfer.h"

// **********
class Ran
{
	public:
	Ran(void)
	{
		srand((unsigned)time(NULL));
	}

// 	Return a random number in the interval [0,1]
	inline double doub()
	{
		return double(rand())/double(RAND_MAX);
	}

//	TODO: Overload doub() function to return a random number in [a,b]

	double doub(double a, double b)
	{
		double r = double(rand())/double(RAND_MAX);
		r *= std::abs(a-b);
		r += MIN(a,b);
		return r;
	}

};

#endif
