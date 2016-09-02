#include <iostream>
#include <cmath>

double GS_recursive(int power)
{
	return (power == 0) ? 1 : (power == 1) ?
	((sqrt(5) - 1) / 2) : (GS_recursive(power - 2) - GS_recursive(power - 1));
}

double GS_exponential(int power)
{
	double phi = (sqrt(5)-1.)/2.;
	double r = 1;
	for (int i = 0; i < power; ++i)
		r *= phi;
	return r;
}

int main ()
{
	std::cout << GS_recursive(50) << std::endl;
	std::cout << GS_exponential(50) << std::endl;
	

}