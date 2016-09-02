#ifndef _THREEDIM_H_
#define _THREEDIM_H_

#include "helfer.h"
#include "adapt.h"

// **********
template <class T>
class Quad3D
{
	struct NRf3
	{
		double xsav, ysav;
		T& func3d;

		NRf3(const Quad3D& q) : func3d(q.func3d)
		{
		}

		double operator() (const double z)
		{
			return func3d(xsav,ysav,z);
		}
	};

	struct NRf2
	{
		T& func3d;
		NRf3 f3;
		const double eps;

		NRf2(const Quad3D& q) : func3d(q.func3d), f3(q), eps(q.eps)
		{
		}

		double operator() (const double y)
		{
			f3.ysav = y;
	
			Adapt a(eps);
			return a.integrate(f3, func3d.z0(f3.xsav,y), func3d.z1(f3.xsav,y));
		}
	};

	struct NRf1
	{
		T& func3d;
		NRf2 f2;
		const double eps;

		NRf1(const Quad3D& q) : func3d(q.func3d), f2(q), eps(q.eps)
		{
		}

		double operator() (const double x)
		{
			f2.f3.xsav = x;

			Adapt a(eps);
			return a.integrate(f2, func3d.y0(x), func3d.y1(x));
		}
	};

// 	Change the required accuracy if necessary
	T& func3d;
	double eps;
	NRf1 f1;
	double x0, x1;

	public:
// 	Call the constructor e.g. as Quad3D<Functor3D> q(func), where func:
// 	1. has public values for the x boundaries as x0, x1;
// 	2. has public functions y0(double) and y1(double) for the y boundaries; and
// 	3. has public functions z0(double,double) and z1(double,double) for the z boundaries
// 	4. must support calls to operator()(double,double,double)

// 	Constructor: get threedimensional integrand and initialize internal objects
	Quad3D(T& func3dIn) : func3d(func3dIn), eps(1e-12), f1(*this)
	{
		x0 = func3d.x0;
		x1 = func3d.x1;
	}

// 	Perform integration for a given accuracy (with a default value 1e-12)
	double integrate(const double epsIn = 1e-12)
	{
		eps = epsIn;
		Adapt a(eps);
		return a.integrate(f1,x0,x1);
	}
};

#endif
