#include <iostream>
#include <cmath>
#include <fstream>

/*double GS_recursive(int power)
{
	return (power == 0) ? 1 : (power == 1) ?
	((sqrt(5) - 1) / 2) : (GS_recursive(power - 2) - GS_recursive(power - 1));
}*/

double GS_recursive(int power)
{
	double* gsvalues = new double[power+1];
	gsvalues[0] = 1.;
	gsvalues[1] = 0.618003398;
	
	for (int i = 2; i <= power; ++i)
	{
		gsvalues[i] = gsvalues[i-2]-gsvalues[i-1];
	}
	
	return gsvalues[power];
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
	std::ofstream file;
	file.open("Golden-Section-plot.txt");
	
	
	for (int i = 1; i <= 50; ++i) 
	{
		file << i << "\t" << GS_recursive(i) << "\t" << i << "\t" << GS_exponential(i) << std::endl;
	}

	
	file.close();
}