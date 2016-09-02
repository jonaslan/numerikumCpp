#ifndef _ITGR_H_
#define _ITGR_H_

#include "helfer.h"

// **********
template <class T>
class IntegrBase
{
	protected:
	T* pde;
	double (T::*f)(const double);

	public:
	IntegrBase(T* pdeIn, double (T::*func)(const double)) : pde(pdeIn), f(func)
	{
	}

	virtual inline double operator() (const double x) = 0;
};

// **********
template <class T>
class Integrand0 : public IntegrBase<T>
{
	using IntegrBase<T>::pde;
	using IntegrBase<T>::f;

	public:
	Integrand0(T* pdeIn, double (T::*func)(const double));
	inline double operator() (const double x);
};

// **********
template <class T>
Integrand0<T>::Integrand0(T* pdeIn, double (T::*func)(const double)) : IntegrBase<T>(pdeIn,func)
{
}

// **********
template <class T>
double Integrand0<T>::operator() (const double x)
{
	return (pde->*f)(x);
}

// **********
template <class T>
class Integrand1 : public IntegrBase<T>
{
	using IntegrBase<T>::pde;
	using IntegrBase<T>::f;
	double c1, c0, d1, d0;

	public:
	Integrand1(T* pdeIn, double (T::*func)(const double), const double c1In, const double c0In, const double d1In, const double d0In);
	inline double operator() (const double x);
};

// **********
template <class T>
Integrand1<T>::Integrand1(T* pdeIn, double (T::*func)(const double), const double c1In, const double c0In, const double d1In, const double d0In) :
	IntegrBase<T>(pdeIn,func)
{
	c1 = c1In;
	c0 = c0In;
	d1 = d1In;
	d0 = d0In;
}

// **********
template <class T>
double Integrand1<T>::operator() (const double x)
{
	return (pde->*f)(x) * ( c1*x + c0 ) * ( d1*x + d0 );
}

#endif
