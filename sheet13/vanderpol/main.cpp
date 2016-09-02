// **************************
// *** Numerikum, WS 2012 ***
// *** Sheet 13 Problem 3 ***
// **************************

#include <iostream>
#include <fstream>
#include "stepperdopr5.h"
#include "vdp.h"
using namespace std;

// **********
int main()
{
// 	Set some necessary parameters
	static const double atol = 1e-14;
	static const double rtol = atol;
	static const double h1 = 0.01, hmin = 0.;

//	Initializing ODE object
	VanDerPol vdp;

//	Initializing storage object for data output
	Output out(1000);

//	Initializing ODE integration object templated on the method
	Odeint<StepperDopr5<VanDerPol>> ode( vdp.ystart, 0.,5., atol,rtol, h1,hmin, out,vdp );

//	Integrate the ODE
	cout << "Starting the integration..." << endl;
	ode.integrate();
	cout << "ODE integration required " << out.steps[0] << " steps." << endl;

//	Write output
	ofstream datei("vdp.txt",ios::trunc);
	for (int i=0; i<out.count; ++i)
		datei << out.xsave[i] << '\t' << out.ysave[0][i] << endl;
	datei.close();
}
