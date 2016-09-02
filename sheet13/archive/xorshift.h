#include <cstdlib>
using namespace std;

#ifdef _MSC_VER
typedef unsigned __int64 Ullong;
#else
typedef unsigned long long int Ullong;
#endif

// **********
class Xorshift
{
	Ullong x;
	double div;

	public:
	Xorshift(Ullong j=rand()) : x(4101842887655102017LL), div(1.)
	{
		x ^= j;
		x = int64();

		div = 1.;
		for (int i=0; i<64; ++i)
			div *= 0.5;
	}

	Ullong int64()
	{
		x ^= x >> 21;
		x ^= x << 35;
		x ^= x >> 4;
		return x * 2685821657736338717LL;
	}

	double doub()
	{
		return div * int64();
	}

	double doub(const double a, const double b)
	{
		return a + (b-a) * doub();
	}
};
