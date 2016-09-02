#ifndef _XORSHIFT
#define _XORSHIFT
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
	Xorshift(Ullong j=rand());
	
	Ullong int64();
	
	double doub();
	
	double doub(const double a, const double b);
};
#endif