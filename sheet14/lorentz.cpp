// Charged particles in a magnetic field
#include <fstream>
#include <iostream>
#include <cmath>
#include "helfer.h"

using namespace std;

// Implementation of Runge-Kutta method of 4th order
class RKSolver
{

	double h; // Step size
	VecDoub k1,k2,k3,k4; // Derivatives of y associated to each of the 4 intermediate points. The full Runge-Kutta step uses a weighted average of them
	

public:
	RKSolver(const double hIn) : h(hIn)
	{
		// Important: the size of the VecDoub objects can only be specified in the constructor, not when initialising the object above (i.e. VecDoub k1(6) in the private definitons would not have worked
		k1.resize(6);
		k2.resize(6);
		k3.resize(6);
		k4.resize(6);
	}

	template <typename T>
	void integrate(T& func,const double t, VecDoub& y)
	{
		VecDoub y1(6),y2(6),y3(6);
		func(t, y, k1); // set the value of k1 
		// Calculating corresponding increment for y, which is used in the evaluation of the next k
		for (int i = 0; i < y.size(); ++i)
			y1[i] = y[i]+0.5*h*k1[i];
		func(t+0.5*h, y1, k2); // set the value of k2

		for (int i = 0; i < y.size(); ++i)
			y2[i] = y[i]+0.5*h*k2[i];
		func(t+0.5*h, y2, k3); // set the value of k3

		for (int i = 0; i < y.size(); ++i)
			y3[i] = y[i]+h*k3[i];
		func(t+h, y3, k4); // set the value of k4			
		
		for (int i = 0; i < y.size(); ++i)
		{
			y[i] = y[i] + h*(1./6.*k1[i]+1./3.*(k2[i]+k3[i])+1./6.*k4[i]);
		}
	}
};

// Class containing the differential equation to be solved
class Lorentz
{
	double B0; // Constant magnetic field
public:
	VecDoub ystart;
	Lorentz() : B0(1.) // As an example, we initalise to 1 the constant magnetic field
	{
		ystart.resize(6); // ystart contains the initial conditions
		// initial position
		ystart[0] = 0.;
		ystart[1] = 0.;
		ystart[2] = 0.;
		
		// initial velocity
		ystart[3] = 1./sqrt(2);
		ystart[4] = 0;
		ystart[5] = 1./sqrt(2);
	}
	// The following overloading operation sets the values of the object dydt according to the RHS of our system of ODE's, i.e. updates them for every iteration
	void operator() (const double t, const VecDoub& y, VecDoub& dydt)
	{
		dydt[0] = y[3];
		dydt[1] = y[4];
		dydt[2] = y[5];
		dydt[3] = y[4]*B0;
		dydt[4] = -y[3]*B0;
		dydt[5] = 0.;
		// Note that i this case we have no dependence neither on t nor in the position (y[0],y[1],y[2]) on the RHS of the equation
	}
};

int main ()
{
	// Open files that will store the resuts of the integration
	ofstream file;
	file.open("lorentz.txt",ios::trunc);
	// Set step size
	double stepsize = 0.1;
	// Instantiate RKSolver object and a functor object for the differential equation
	RKSolver RKS(stepsize);
	Lorentz L;
	 
	double t = 0.; // Time variable
	const double tmax = 1000000.; // Upper limit of integration
	VecDoub y = L.ystart; // Initialise the object with the initial conditions
	// Initial velocity calculation: since there's only a magnetic field, the kinetic energy of the particle remains constant. The divergence of the magnitude of the velocity, v, from the initial velocity, v0, as the integration proceeds is a measure of the convergence of the numerical method
	double v0 = sqrt(y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
	double v = v0;
	
	file << "#t\tx\ty\tz\t|v|" << endl;
	
	while (t < tmax && (((v/v0) < 1.01)&&((v/v0) > 0.99)))
	{
		// Write current position and velocity to file
 
		file << t << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << v << endl;
		RKS.integrate<Lorentz>(L,t,y); // Integration step
		t = t + stepsize; // Update time
		v = sqrt(y[3]*y[3]+y[4]*y[4]+y[5]*y[5]); // Calculation of magnitude of velocity
	}
	cout << t << "\t" << v << endl;	
}
