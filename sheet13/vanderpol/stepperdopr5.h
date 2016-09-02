// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 13 Problem 3 ***
// **************************

#ifndef _STEPPERDOPR5_H_
#define _STEPPERDOPR5_H_

#include "helfer.h"
#include "odeint.h"
#include "vdp.h"

// **********
template <class D>
struct StepperDopr5 : StepperBase
{
	typedef D Dtype;
	VecDoub k2,k3,k4,k5,k6;
	VecDoub rcont1,rcont2,rcont3,rcont4,rcont5;
	VecDoub dydxnew;

	StepperDopr5(VecDoub& yy, VecDoub& dydxx, double& xx, const double atoll, const double rtoll, bool dens);

	void step(const double htry, D& derivs);
	void dy(const double h, D& derivs);
	void prepare_dense(const double h, D& derivs);
	double dense_out(const int i, const double x, const double h);
	double error();

	struct Controller
	{
		double hnext, errold;
		bool reject;
		Controller();
		bool success(const double err, double& h);
	};

	Controller con;
};

#endif
