#include "xorshift.h"

Xorshift::Xorshift(Ullong j) : x(4101842887655102017LL), div(1.)
{
	x ^= j;
	x = int64();

	div = 1.;
	for (int i=0; i<64; ++i)
		div *= 0.5;
}

Ullong Xorshift::int64()
{
	x ^= x >> 21;
	x ^= x << 35;
	x ^= x >> 4;
	return x * 2685821657736338717LL;
}

double Xorshift::doub()
{
	return div * int64();
}

double Xorshift::doub(const double a, const double b)
{
	return a + (b-a) * doub();
}
