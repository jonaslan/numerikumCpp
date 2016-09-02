#ifndef _INTERP_H_
#define _INTERP_H_

#include "helfer.h"
using namespace std;

// **********
class Base_interp
{
	protected:
	int n, mm, jsav, cor, dj;
	const double *xx, *yy;

	public:
	Base_interp(VecDoub& x, const double *y, int m)
		: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
	{
		dj = MIN(1,(int)pow((double)n,0.25));
	}

	double interp(double x)
	{
		int jlo = cor ? hunt(x) : locate(x);
		return rawinterp(jlo,x);
	}

	int locate(const double x);
	int hunt(const double x);
	
	virtual double rawinterp(int jlo, double x) = 0;
};

// **********
int Base_interp::locate(const double x)
{
	int ju,jm,jl;

	if (n < 2 || mm < 2 || mm > n)
		throw("locate size error");

	bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;

	while (ju-jl > 1)
	{
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}

	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;

	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

// **********
int Base_interp::hunt(const double x)
{
	int jl=jsav, jm, ju, inc=1;

	if (n < 2 || mm < 2 || mm > n)
		throw("hunt size error");

	bool ascnd=(xx[n-1] >= xx[0]);

	if (jl < 0 || jl > n-1)
	{
		jl=0;
		ju=n-1;
	}
	else
	{
		if (x >= xx[jl] == ascnd)
		{
			while (true)
			{
				ju = jl + inc;

				if (ju >= n-1)
				{
					ju = n-1;
					break;
				}

				else if (x < xx[ju] == ascnd)
					break;

				else
				{
					jl = ju;
					inc += inc;
				}
			}
		}
		else
		{
			ju = jl;
			while (true)
			{
				jl = jl - inc;

				if (jl <= 0) {
					jl = 0;
					break;
				}

				else if (x >= xx[jl] == ascnd)
					break;

				else
				{
					ju = jl;
					inc += inc;
				}
			}
		}
	}

	while (ju-jl > 1)
	{
		jm = (ju+jl) >> 1;

		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}

	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;

	return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

// **********
class Poly_interp : public Base_interp
{
	double dy;

	public:
	Poly_interp(VecDoub& xv, VecDoub& yv, int m) : Base_interp(xv,&yv[0],m), dy(0.)
	{
	}

	double rawinterp(int jl, double x);
};

// **********
double Poly_interp::rawinterp(int jl, double x)
{
	int i,m,ns=0;
	double y,den,dif,dift,ho,hp,w;
	const double *xa = &xx[jl], *ya = &yy[jl];

	VecDoub c(mm),d(mm);
	dif=abs(x-xa[0]);

	for (i=0;i<mm;i++)
	{
		if ((dift=abs(x-xa[i])) < dif)
		{
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	y=ya[ns--];

	for (m=1;m<mm;m++)
	{
		for (i=0;i<mm-m;i++)
		{
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];

			if ((den=ho-hp) == 0.0)
				throw("Poly_interp error");

			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	}
	return y;
}

#endif
