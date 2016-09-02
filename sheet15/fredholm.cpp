// Fredholm

#include <iostream>
#include <cmath>
#include "helfer.h"
#include "gauleg.h"
#include "ludcmp.h"


using namespace std;


class Fredh2
{
	int N; //Number of points in which we discretise the integration from a to b
	double a,b,lambda; // lambda: parameter of the equation
	VecDoub fvec,t,w;
	
public:
	
	Fredh2(int n, double lambdaIn) : N(n), a(0.), b(M_PI), lambda(lambdaIn)
	{
		fvec.resize(N);
		t.resize(N);
		w.resize(N); // t, fvec and w are in principle empty vectors

		// Determine the weights corresponding to the different points in which we have discretised our integral using Gauss-Legendre method
		// Shift interval to: a to b(remember Legendre polynomials are defined between -1 and 1)
		gauleg(a,b,t,w);

		// Declaration of matrix mat and vector gvec needed to calculate the values of the function for the determined abcissas
		MatDoub mat(N,N);
		VecDoub gvec(N);
		
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				// Build matrix mat
				if (i == j)
					mat[i][j] = 1.-lambda*w[j]*(t[i]-t[j]);
				else
					mat[i][j] = -lambda*w[j]*(t[i]-t[j]);
			}
			// Build vector gvec 
			gvec[i] = sin(t[i]);
		}
		// Instantiate object to solve matrix equation
		LUdcmp alu(mat);
		// Solve matrix equation and determine fvec
		alu.solve(gvec,fvec);
	}
	// This overloaded operator interpolates the solution of the integral equation for arbitrary t's
	double operator () (const double tIn)
	{
		double sum = 0.;
		for (int j = 0; j < N; ++j)
		{
			sum += w[j]*(tIn-t[j])*fvec[j];
		}
		return lambda*sum+sin(tIn);
	}
};


int main ()
{
	
	double lambda = 0.1;
	int N = 100; // Integration points
	Fredh2 Fred(N,lambda); // Integral equation solver instantiation

	// Initialisation of constants for analytical calculation
	double A1, A2,g1,g2,delta1,delta2,delta3;
	g1 = 2.;
	g2 = M_PI;
	delta1 = M_PI;
	delta2 = M_PI*M_PI;
	delta3 = M_PI*M_PI*M_PI;
	A1 = (12.*g1 + 6.*lambda*(g1*delta2-2.*g2*delta1))/(lambda*lambda*delta1*delta1*delta1*delta1+12.);
	A2 = (-12.*g2 + 2.*lambda*(3.*g2*delta2-2.*g1*delta3))/(lambda*lambda*delta1*delta1*delta1*delta1+12.);

	// Variables to check our solution
	int n = 100; // Number of points on which to evaluate the solved f(t) (for both analytical and numerical solution)
	double t_0 = 0; // Start of our plot
	double t_f = M_PI; // End of our plot (note that this t_0 and t_f can in principle be completely different from a and b (limits of integration in the integral equation)
	double step = (t_f - t_0)/(n-1);
	double t = t_0;

	cout << "# t f(t)_Analytical f(t)_numerical" << endl;
	cout << "#t_0 = " << t_0 << "\tt_f = " << t_f << endl;

	// Gather results of analytical and numerical calculations
	for (int i = 0; i < n; ++i) 
	{
		double f = sin(t)+lambda*(A1*t+A2); // Analytical solution

		// Cout-ing results to standard output, which can be saved to a file through ">" when calling the program
		cout << t << "\t" << f << "\t" << Fred(t) << endl;
		t += step;
	}	
}
