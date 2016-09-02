#include <iostream>
#include <string>
#include <limits>
#include <sstream>
#include <cmath>


std::string checkNum(double i)
{
	if (i != i)
		return "not a number\n";
	if (!std::isfinite(i))
		return "infinite\n";
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
	double pInf = log(+0);
	double nInf = -1/.0;
	std::cout << checkNum(d) << std::endl;
	std::cout << checkNum(pInf) << std::endl;
	std::cout << checkNum(nInf) << std::endl;
	
	if (pInf==nInf)
	{
		std::cout << "Negative inf = Positive inf" << std::endl;
	} else {
		std::cout << "Negative inf != Positive inf" << std::endl;
	}

	if (d==pInf)
	{
		std::cout << "Not A Number = Positive inf" << std::endl;
	} else {
		std::cout << "Not A Number != Positive inf" << std::endl;
	}
	
	std::cout << std::cout.precision(20) << epsilon() << std::endl;
}