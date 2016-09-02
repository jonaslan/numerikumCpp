#include <iostream>
#include <limits>
#include <cmath>
#include <fstream>


// Windows/Linux-compatible definition of unsigned long long int
#ifdef _MSC_VER
typedef unsigned __int64 Ullong;
#else
typedef unsigned long long int Ullong;
#endif



class Random
{
	// private variables
	Ullong seed,v,a,m;
	double ran;
	int a1,a2,a3,c;
	
	// private xorshift function.This will never be directly called by the users
	// and can therefore remain private. It takes unsigned long long, shifts it and xors it
	// using three private variables and returns input to the randval MLCG
	Ullong xorshift(Ullong x)
	{
		x = x^(x >> a1);
		x = x^(x << a2);
		x = x^(x >> a3);
		return x;
	}

public:
	// constructor initializes all private variables and takes input from user which
	// is used to xor out the seed for the xorshift function
	Random(const int& input)
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
	
	// the actual output function called by the user. A multiplicative linear
	// congruential generator which gets its input from the xorshift. The seed is set
	// to the output.
	Ullong randval()
	{
		seed = (a*xorshift(seed)+c)%m;
		return seed;
	}
	
	// the inrange function returns a random value in a specific range.
	// It first casts the output from the randval function and the maximum value for
	// the unsigned long long as doubles and normalizes the random value to a value
	// between 0 and 1. Then it multiplies with the range and adds the lower of the
	// range values
	double inrange(const double& a, const double& b)
	{	
		ran = static_cast<double>(randval())/static_cast<double>(std::numeric_limits<Ullong>::max());
		ran *= std::abs(a-b);
		return (a>b)?(b+ran):(a+ran);
	}

};

int main ()
{
	
	// Tests and user input
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

	// File outputs for the plots
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