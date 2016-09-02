// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 13 Problem 3 ***
// **************************

#ifndef DEINT_H_
#define DEINT_H_

#include <climits>
#include "helfer.h"

// **********
class StepperBase
{
	public:
	double &x;
	double xold;
	VecDoub &y, &dydx;
	double atol, rtol;
	bool dense;
	double hdid;
	double hnext;
	double EPS;
	int n, neqn;
	VecDoub yout, yerr;

// 	Parent class for stepper classes StepperBS and StepperDopr853
	StepperBase(VecDoub& yy, VecDoub& dydxx, double& xx, const double atoll, const double rtoll, bool dens) :
		x(xx), y(yy), dydx(dydxx), atol(atoll), rtol(rtoll), dense(dens), n(y.size()), neqn(n), yout(n), yerr(n)
	{
	}
};

// **********
struct Output
{
	int kmax;
	int nvar;
	int nsave;
	bool dense;
	int count;
	double x1,x2,xout,dxout;

	inline void init(const int neqn, const double xlo, const double xhi);
	inline void resize(void);

	template <class Stepper>
	inline void save_dense(Stepper& s, const double xout, const double h);
	inline void save(const double x, VecDoub& y);

	template <class Stepper>
	inline void out(const int nstp, const double x, VecDoub& y, Stepper& s, const double h);

	public:
	int steps[2];

	Output(void);
	Output(const int nsavee);

	VecDoub xsave;
	MatDoub ysave;
};

// **********
template <class Stepper>
class Odeint
{
	friend class Output;

	static const int MAXSTP=1000000000;
	double EPS;
	int nvar;
	double x1,x2,hmin;
	bool dense;
	VecDoub y,dydx;
	VecDoub &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	int nstp;
	double x,h;

	public:
	int nok;
	int nbad;

	// Constructor
	Odeint(VecDoub& ystartt, const double xx1, const double xx2, const double atol, const double rtol, \
	const double h1, const double hminn, Output& outt, typename Stepper::Dtype& derivss);

	// Integration routine
	void integrate(void);
};

#endif
