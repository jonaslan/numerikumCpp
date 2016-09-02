#include<cmath>
#include "helfer.h"
#include<iostream>
#include<fstream>

// a power function without the heavy bells and whistles of pow
static inline double power(const double& x, const int& p)
{
	double r = x;
	for (int i = 1; i < p; ++i)
		r *= x;
	return r;
}



class Dawson
{
private:
	double factor1,factor2,factor3,h,n0,N,x0,xp,invpi;

    // Private Taylor Series Expansion function
	double Taylor(const double& x)
	{
	    return x - (2./3.)*power(x,3) + (4./15.)*power(x,5) - (8./105.)*power(x,7);
	}

public:
    // Constructor with initialization list for all members which don't require input data
	Dawson() : h(0.4), N(11), invpi(1./sqrt(M_PI)), factor2(exp(-h*h)) {}

	double calculate(const double& xin)
	{
        // Sign change if necessary.
	    double x = (xin>0) ? xin : -xin;

	    // If x is too small, call Taylor Series approximation
	    if (x < 0.2)
        {
            return Taylor(x);
        }

        // Define n0. If the integer above is divisible by 2, n0 assumes its value. Else, the value below.
        // Example: x = 3.4. Function checks 4%2 == 0 and evaluates to true -> n0 = 4
        //          x = 2.9. Function checks 3%2==0 and evaluates to false -> n0 -> 2
        if (!((int)ceil(x/h)%2))
        {
            n0 = ceil(x/h);
        }
        else
            n0 = floor(x/h);

        // Assignment of internal factors
        x0 = n0*h;
        xp = x - x0;
        factor1 = exp(-power(xp,2));
        factor3 = exp(2.*xp*h);

        // Summation of all terms in the span
        double sum = 0;
        for (int n = -N; n <= N; n = n + 2)
        {
            sum += (1./(n+n0))*power(factor2,n*n)*power(factor3,n);
        }

        // Return sum with correct sign.
        if (xin < 0)
            return -invpi*factor1*sum;
        else
            return invpi*factor1*sum;
	}
};

int main ()
{
    double upper = 0.54;
    double lower = -0.54;
	double x = lower;
	Dawson D;
    double ds[100];
    double xs[100];

    std::ofstream file;
    file.open("Dawson-values.txt",std::ios::trunc);

	for (int i = 0; i < 1000; ++i)
    {
        file << x << "\t" << D.calculate(x) << std::endl;
        std::cout << D.calculate(x) << std::endl;
        x += (upper-lower)/1000.;
    }

    file.close();
}
