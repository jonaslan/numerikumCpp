#include <iostream>
#include <string>
#include <limits>
#include <sstream>
#include <cmath>


std::string checkNum(double i)
{
	if (i != i)
		return " is not a number\n";
	if (!std::isfinite(i))
		return " is infinite\n";
	return 0;
}

double epsilon()
{
	volatile double diff = 0.;
	int count = 0;
	double eps = 2e-60;
	
	while ( diff == 0 && count <= 1000)
	{
		++count;
		diff = (1.+eps) - 1.;
		eps *= 2;
	}
	return eps;
}


int main ()
{
	double d = sqrt(-1);
	double e = log(+0);
	std::cout << d << checkNum(d) << std::endl;
	std::cout << e << checkNum(e) << std::endl;
	std::cout << std::cout.precision(20) << epsilon() << std::endl;
}