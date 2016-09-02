#include <iostream>
#include "fem1d.h"
#include "ludcmp.h"
#include "adapt.h"
#include "integrand.h"

using namespace std;

// Class containing the methods a,b,c and f corresponding to the coefficients of the solver implemented in the header "fem1d.h"
class Function
{
	double s; // exponent that determines the kind of turbulence (0 for isotropic, 2/3 for 1D)

public:
	
	double m0; // initial value for pitch-angle cosine

	Function() : s(0.), m0(0.7)
	{}
	
	//Implementation of the methods corresponding to each coefficient (note that c and f are exchanged w.r.t. notation in exercise sheet)
	double a (const double x)
	{
		return (1.-x*x)*pow(abs(x),s);
	}

	double b (const double x)
	{
		return 0.;
	}
	
	double c (const double x)
	{
		return 1.;
	}
	
	double f (const double x)
	{
		return 0.;
	}
};


int main ()
{
	Function func;
	
	int Nx = 200.; // Angular resolution (number of points for variable mu, e.g. used to determine f(mu, t = t_0) at each time step) 
	int Mt = 100.;	// Time resolution (number of time steps)
	double dt = 1e-2; // Time increment
	// Declaration of vectors needed for the method Fem1D::run (see header fem1d.h for explanation)
	VecDoub norm, meansq;

	// Instantiate the solver object
	Fem1D<Function> fem(func,Nx,Mt,dt);

	// Start solver
	fem.run(norm,meansq);

/* We used the following code to check the deviation from 1 of the integral of f(mu,t) over mu at each time step t. In all the cases analysed the value of the integral was 1, which is a good sign of the stability of the method 
	for (int i = 0; i < norm.size(); ++i)
	{
		cout <<  "Norm: " << norm[i] << "\tMean squares: " << meansq[i] << endl;
	}
*/	
	cout << "Basic qualitative difference between isotropic and 1D-turbulence:" << endl;
	cout << "In the isotropic case, after a bunch of particles arrives in a certain direction (not necessarily aligned with the magnetic field, which corresponds to mu=1), the isotropic turbulence will induce a diffusive process so that particles will tend to align their \"axis of rotations\" with the homogeneous magnetic field (maximum of the distribution after a certain time is at mu=1, i.e. in the direction of the field). As time passes, however, there will be an increasing portion of particles that will even scatter in directions opposite to the field, and the distribution tends to a constant value over the whole range of mu's for longer times (diffusion in all directions)." << endl;
	cout << "On the contrary, when there is 1d-turbulence in the system, it is directed in a direction parallel to the field, so that it creates a kind of \"potential barrier\" such that particles cannot easily change their initial directions. This is seen in the two plots we show, for opposite initial directions of our beam of particles. It is possible to observe a \"bump\" in mu=0 after the diffusion process occurs for a certain time, which is a sign of the rapid decrease of the distribution when considering components of the parallel velocity which have changed (diffused) and are now opposed to the initial parallel velocity of the particle (unfavourable)" << endl;
}
