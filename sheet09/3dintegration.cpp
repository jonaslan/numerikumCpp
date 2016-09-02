#include <cmath>
#include <iostream>
#include <iomanip> // std::setprecision
#include "threedim.h"

using namespace std;

class Functor3D
{
	public:
		double a, b, c; // Parameters that determine shape of the ellipsoid and the cone
		double x0, x1; // Integration limits for the outer integral (the one that is integrated the last)

		// Constructor. Paramters a, b, c initialised according to problem sheet
		Functor3D() : a(1.), b(2.), c(3.)
	{
		x0 = - a * sqrt(0.5); // These values come from analytical manipulation
		x1 = a * sqrt(0.5); 
	}
		// The following four functions represent the integration limits for the inner variables of integration y and z
		double y0(const double x)
		{
			if ((x/a) * (x/a) >= 0.5) 
			{
				return 0; /* This condition is necessary whenever the argument of sqrt can be negative. This can only happen (when analytical manipulation is done correctly) in the evaluation corresponding to the limits x0 and x1 due to round-off errors */
			}
			else
			{
				return -b * sqrt(0.5-(x/a)*(x/a)); //from analytical manipulation
			}

		}

		double y1(const double x)
		{
			if ((x/a) * (x/a) >= 0.5) 
			{
				return 0;
			}
			else
			{
				return b * sqrt(0.5-(x/a)*(x/a)); //from analytical manipulation
			}

		}


		double z0(const double x, const double y)
		{
			return c * sqrt((x/a)*(x/a) + (y/b) * (y/b)); //from analytical manipulation

		}

		double z1(const double x, const double y)
		{
			if ((x/a) * (x/a) + (y/b) * (y/b) >= 1)
			{
				return 0;
			}
			else
			{
				return c * sqrt(1 - (x/a) * (x/a) - (y/b) * (y/b)); //from analytical manipulation
			}

		}

		// Overloading operator is the three-dimensional integrand. In this case, 1 because we want to calculate only the integral over dV (which gives the total volume)
		double operator () (double x, double y, double z)
		{
			return 1.;
		}


};


int main ()
{
	// Instantiation of Functor3D and Quad3D (function containing methods to calculate volume integral) objects
	Functor3D func;
	Quad3D<Functor3D> q(func);
	// Calculation of the analytical (exact) and numerical value of the integral
	double analytical = M_PI * func.a * func.b * func.c * (2. - sqrt(2))/3.;
	double numerical = q.integrate();
	//Comparison of numerical vs analytical values
	cout << "Numerical approximation to the integral: " << setprecision (15) << numerical << endl;
	cout << "Analytical value of the integral: " << setprecision (15) << analytical << endl;
	/* Numerical approximation to the integral: 3.68060473804278
	   Analytical value of the integral: 3.68060473804244
	   Therefore, the difference between the two values is of the order of 10e-12  */

}

