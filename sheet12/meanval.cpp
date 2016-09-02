#include "xorshift.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include "helfer.h"

// Class that calculates different statistical parameters of a sample
class Distribution
{
private:
	VecDoub values; //Data points
	int N; //Number of data points
	
public:
	Distribution(VecDoub v) : values(v), N( values.size() ) {}
	
	double mean();
	double variance();
	double stddev();
	double skewness();
	double kurtosis();
};

double Distribution::mean()
{
	double m = 0;
	for (int i = 0; i < N; ++i)
		m += values[i];
	return m/N;
}

double Distribution::variance()
{
	double m = mean();
	double v1 = 0;
	double v2 = 0;
	
	for (int i = 0; i < N; ++i)
	{
		v1 += (values[i]-m)*(values[i]-m);
		v2 += (values[i]-m);
	}
	return 1./(N-1)*(v1-1./N*v2*v2); // Note N-1 in the denominator since is the variance obtained from a sample
}

double Distribution::stddev() 
{
	return sqrt(variance());
}

double Distribution::skewness()
{
	double s = 0;
	double m = mean();
	double std = stddev();
	
	for (int i = 0; i < N; ++i)
		s += ((values[i]-m)/std)*((values[i]-m)/std)*((values[i]-m)/std);
	return s/N;
}

double Distribution::kurtosis()
{
	double s = 0;
	double m = mean();
	double std = stddev();
	
	for (int i = 0; i < N; ++i)
		s += ((values[i]-m)/std)*((values[i]-m)/std)*((values[i]-m)/std)*((values[i]-m)/std);
	return s/N-3;
}

// The following function assigns each value to the correspoding bin
int bin_number(const double x, const double max, const double min)
{
	const int M = 256; // Total number of bins 
	int num = 0; // Bin counter
	double v = min;
	
	const double delta = (max-min)/M; //Width of each bin
	
	while(v+delta < x)
	{
		v += delta; 
		++num;
	}
	return num;
}
	

int main()
{
	// Initialisation of random number generator
	std::srand((unsigned) time(NULL));
	Xorshift x;
	
	
	const int N = 1e6; // Number of generated points
	const int M = 256; // number of bins
	

	// Declaration of vectors storing the generated values and 
	VecDoub dist(N,0.);
	VecDoub bins(M,0.);

	
	int count = 0; // Counter for generated points
	
	double eps = numeric_limits<double>::epsilon();	
	
	double max = 0.;
	double min = 0.;
	
	// Generation of random numbers according to algorithm in the instructions. Two random numbers are generated at a time, x1 and x2, and stored in dist
	while (count < N)
	{
		double u1 = x.doub(-1.,1.);
		double u2 = x.doub(-1.,1.);	
		double q = u1*u1+u2*u2;
		if (q <= eps || q >= 1-eps)
			continue;
		
		double p = sqrt(-2*log(q)/q);
			
		double x1 = p*u1;
		double x2 = p*u2;
		
		// Calculation of maximum and minimum of the set of random numbers generated so far
		if (x1 > max)
			max = x1;
		else if (x1 < min)
			min = x1;
		else if (x2 > max)
			max = x2;
		else if (x2 < min)
			min = x2;
		
		
		dist[count++] = x1;
		dist[count++] = x2;
	}
	// Assign each random number generated to the corresponding bin
	for (int i = 0; i < N; ++i)
		++bins[bin_number(dist[i],max,min)];

	// Declaration of a vector that will contain the midpoints of each bin (used later for plotting)
	VecDoub midpoints(M);
 	double delta = (max-min)/M;
	
	midpoints[0] = min+delta/2.;

	// Calculation of distribution parameters
	Distribution stats(dist);

	// Assign midpoints to each bin
	for (int i = 1; i < midpoints.size(); ++i)
		midpoints[i] = midpoints[i-1]+delta;
	cout << "#\t Mean: " << stats.mean() << 
	"\n#\tVariance: " << stats.variance()<< 
	"\n#\tStandard deviation: " << stats.stddev() <<
	"\n#\tSkewness: " << stats.skewness() <<
	"\n#\tKurtosis: " << stats.kurtosis() << std::endl;
	
	// Display counts for each bin 
	for (int i = 0; i < 256; ++i)
		std::cout << midpoints[i] << "\t" << bins[i] << std::endl;

	// Conclusion of analysis: since skewness and kurtosis are very close to 0, we have very good reasons to confirm that the random distribution is Gaussian as expected (mean and variance are also very close to 0 and 1, respectively, which we guess they are the parameters predicted for such a random number generation method)
}
