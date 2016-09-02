#include <iostream>
#include <limits>
#include <cstdlib>
#include <fstream>

#ifdef _MSC_VER
typedef unsigned __int64 Ullong;
#else
typedef unsigned long long int Ullong;
#endif

class Random
{
Ullong seed,v,a,m;
double ran;
int a1,a2,a3,c;

	Ullong xorshift(Ullong x)
	{
		x = x^(x >> a1);
		x = x^(x << a2);
		x = x^(x >> a3);
		return x;
	}

public:
	Random(int input)
	{
		v = 4101842887655102017LL;
		a = 2685821657736338717LL;
		m = 2e64;
		c = 0;
		a1 = 21;
		a2 = 35;
		a3 = 4;
		seed = input^v;
	}
	Ullong randval()
	{
		seed = (a*xorshift(seed)+c)%m;
		return seed;
	}
	
	double inrange(double a, double b)
	{	
		ran = static_cast<double>(randval())/static_cast<double>(std::numeric_limits<Ullong>::max());
		ran *= abs(a-b);
		return (a>b)?(b+ran):(a+ran);
	}

};

int main ()
{
	int s;
	std::cout << "Enter seed value: ";
	std::cin >> s;
	Random R(s);
	double low, high;
	std::cout << "Specify range: \nLower Bound: ";
	std::cin >> low;
	std::cout << "Upper bound: ";
	std::cin >> high;
	std::cout << R.inrange(high,low) << std::endl;

	std::ofstream file;
	file.open("CRNG-Spectral-Test-2D.txt",std::ios::trunc);
	

	for (int i = 0; i < 30000; ++i)
		file << R.randval() << "\t" << R.randval() << std::endl;
	
	file.close();
	file.open("CRNG-Spectral-Test-3D.txt",std::ios::trunc);
	
	for (int i = 0; i < 30000; ++i)
		file << R.randval() << "\t" << R.randval() << "\t" << R.randval() << std::endl;
		
	file.close();
}