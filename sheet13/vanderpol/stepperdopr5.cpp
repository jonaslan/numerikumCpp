// **************************
// *** Numerikum, WS 2013 ***
// *** Sheet 13 Problem 3 ***
// **************************

#include "stepperdopr5.h"

// **********
template <class D>
StepperDopr5<D>::StepperDopr5(VecDoub& yy, VecDoub& dydxx, double& xx, const double atoll, const double rtoll, bool dens) :
	StepperBase(yy,dydxx,xx,atoll,rtoll,dens), k2(n), k3(n), k4(n), k5(n), k6(n), \
	rcont1(n), rcont2(n), rcont3(n), rcont4(n), rcont5(n), dydxnew(n)
{
	EPS = numeric_limits<double>::epsilon();
}

// **********
template <class D>
void StepperDopr5<D>::step(const double htry, D& derivs)
{
	double h=htry;

	for (;;)
	{
		dy(h,derivs);
		double err = error();
		if (con.success(err,h))
			break;
		if (fabs(h) <= fabs(x)*EPS)
			throw("Stepsize underflow in StepperDopr5");
	}

	if (dense)
		prepare_dense(h,derivs);

	dydx = dydxnew;
	y = yout;
	xold = x;
	x += (hdid=h);
	hnext = con.hnext;
}

// **********
template <class D>
void StepperDopr5<D>::dy(const double h, D& derivs)
{
	static const double c2=0.2, c3=0.3, c4=0.8, c5=8./9., a21=0.2, a31=3./40.;
	static const double a32=9./40., a41=44./45., a42=-56./15., a43=32./9., a51=19372./6561.;
	static const double a52=-25360./2187., a53=64448./6561., a54=-212./729., a61=9017./3168.;
	static const double a62=-355./33., a63=46732./5247., a64=49./176., a65=-5103./18656.;
	static const double a71=35./384., a73=500./1113., a74=125./192., a75=-2187./6784.;
	static const double a76=11./84., e1=71./57600., e3=-71./16695., e4=71./1920.;
	static const double e5=-17253./339200., e6=22./525., e7=-1./40.;

	VecDoub ytemp(n);

	for (int i=0; i<n; i++)
		ytemp[i] = y[i]+h*a21*dydx[i];
	derivs(x+c2*h,ytemp,k2);

	for (int i=0; i<n; i++)
		ytemp[i] = y[i]+h*(a31*dydx[i]+a32*k2[i]);
	derivs(x+c3*h,ytemp,k3);

	for (int i=0; i<n; i++)
		ytemp[i] = y[i]+h*(a41*dydx[i]+a42*k2[i]+a43*k3[i]);
	derivs(x+c4*h,ytemp,k4);

	for (int i=0; i<n; i++)
		ytemp[i] = y[i]+h*(a51*dydx[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
	derivs(x+c5*h,ytemp,k5);

	for (int i=0; i<n; i++)
		ytemp[i] = y[i]+h*(a61*dydx[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
	double xph = x+h;
	derivs(xph,ytemp,k6);

	for (int i=0; i<n; i++)
		yout[i] = y[i]+h*(a71*dydx[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
	derivs(xph,yout,dydxnew);

	for (int i=0; i<n; i++)
		yerr[i] = h*(e1*dydx[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydxnew[i]);
}

// **********
template <class D>
void StepperDopr5<D>::prepare_dense(const double h, D& derivs)
{
	VecDoub ytemp(n);
	static const double d1=-12715105075./11282082432., d3=87487479700./32700410799.;
	static const double d4=-10690763975./1880347072., d5=701980252875./199316789632.;
	static const double d6=-1453857185./822651844., d7=69997945./29380423.;

	for (int i=0;i<n;i++)
	{
		rcont1[i] = y[i];
		double ydiff = yout[i]-y[i];
		rcont2[i] = ydiff;
		double bspl = h*dydx[i]-ydiff;
		rcont3[i] = bspl;
		rcont4[i] = ydiff-h*dydxnew[i]-bspl;
		rcont5[i] = h*(d1*dydx[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*dydxnew[i]);
	}
}

// **********
template <class D>
double StepperDopr5<D>::dense_out(const int i, const double x, const double h)
{
	double s=(x-xold)/h;
	double s1=1.-s;
	return rcont1[i] + s*( rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*rcont5[i])) );
}

// **********
template <class D>
double StepperDopr5<D>::error()
{
	double err=0.0, sk;
	for (int i=0; i<n; i++)
	{
		sk = atol + rtol * MAX( fabs(y[i]),fabs(yout[i]) );
		err += SQR(yerr[i]/sk);
	}
	return sqrt(err/n);
}

// **********
template <class D>
StepperDopr5<D>::Controller::Controller() : reject(false), errold(1e-4)
{
}

// **********
template <class D>
bool StepperDopr5<D>::Controller::success(const double err, double& h)
{
	static const double beta=0., alpha=0.2-beta*0.75, safe=0.9, minscale=0.2, maxscale=10.;
	double scale;

	if (err <= 1.)
	{
		if (err == 0.)
			scale = maxscale;
		else
		{
			scale = safe * pow(err,-alpha)*pow(errold,beta);
			if (scale < minscale)
				scale = minscale;
			if (scale > maxscale)
				scale = maxscale;
		}
		if (reject)
			hnext = h*scale;

		errold = MAX(err,1.0e-4);
		reject = false;
		return true;
	}

	else
	{
		scale = MAX(safe*pow(err,-alpha),minscale);
		h *= scale;
		reject = true;
		return false;
	}
}

template class StepperDopr5<VanDerPol>;
