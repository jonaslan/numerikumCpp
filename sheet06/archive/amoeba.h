#ifndef _AMOEBA_H_
#define _AMOEBA_H_

#include "helfer.h"

// **********
class Amoeba
{
	const double ftol;
	int nfunc, mpts, ndim;

	VecDoub y;
	MatDoub p;

	public:
// 	Constructor: initialize Amoeba object with required accuracy
	Amoeba(const double ftoll = 1e-15) : ftol(ftoll)
	{
	}

// 	Minimization function to be called by user:
//	return value: vector containing the coordinates of the minimum
// 	point = Vector object containing the starting value
// 	del = initial step size (use, e.g., 0.1)
// 	fmin = function value at minimum will be stored here
// 	func = function/functor to be minimized
	template <class T>
	VecDoub minimize(VecDoub& point, const double del, double& fmin, T& func)
	{
		VecDoub dels(point.size(),del);
		return minimize(point, dels, fmin, func);
	}

	private:
	template <class T>
	VecDoub minimize(VecDoub& point, VecDoub& dels, double& fmin, T& func)
	{
		int ndim=point.size();
		MatDoub pp(ndim+1,ndim);

		for (int i=0;i<ndim+1;i++)
		{
			for (int j=0;j<ndim;j++)
				pp[i][j]=point[j];
			if (i !=0 )
				pp[i][i-1] += dels[i-1];
		}
		return minimize(pp, fmin, func);
	}

	template <class T>
	VecDoub minimize(MatDoub& pp, double& fmin, T& func)
	{
		const int NMAX=5000;
		const double TINY=1.0e-10;
		int ihi,ilo,inhi;
		mpts = pp.nrows();
		ndim = pp.ncols();
		VecDoub psum(ndim), pmin(ndim), x(ndim);
		p=pp;
		y.resize(mpts);

		for (int i=0;i<mpts;i++)
		{
			for (int j=0;j<ndim;j++)
				x[j]=p[i][j];
			y[i]=func(x);
		}

		nfunc=0;
		get_psum(p,psum);

		int iter = 0;
		while (true)
		{
			++iter;
			ilo=0;
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);

			for (int i=0;i<mpts;i++)
			{
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi])
				{
					inhi=ihi;
					ihi=i;
				}
				else if (y[i] > y[inhi] && i != ihi)
					inhi=i;
			}

			double rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);

			if (rtol < ftol)
			{
				SWAP(y[0],y[ilo]);
				for (int i=0;i<ndim;i++)
				{
					SWAP(p[0][i],p[ilo][i]);
					pmin[i]=p[0][i];
				}
				fmin=y[0];
				return pmin;
			}

			if (nfunc >= NMAX)
				throw("NMAX exceeded");

			nfunc += 2;
			double ytry=amotry(p,y,psum,ihi,-1.0,func);

			if (ytry <= y[ilo])
				ytry=amotry(p,y,psum,ihi,2.0,func);

			else if (ytry >= y[inhi])
			{
				double ysave=y[ihi];
				ytry=amotry(p,y,psum,ihi,0.5,func);

				if (ytry >= ysave)
				{
					for (int i=0;i<mpts;i++)
					{
						if (i != ilo)
						{
							for (int j=0;j<ndim;j++)
								p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
							y[i]=func(psum);
						}
					}
					nfunc += ndim;
					get_psum(p,psum);
				}
			}

			else
				--nfunc;
		}
	}

	inline void get_psum(MatDoub& p, VecDoub& psum)
	{
		for (int j=0;j<ndim;j++)
		{
			double sum=0.0;
			for (int i=0;i<mpts;i++)
				sum += p[i][j];
			psum[j]=sum;
		}
	}

	template <class T>
	double amotry(MatDoub& p, VecDoub& y, VecDoub& psum, const int ihi, const double fac, T& func)
	{
		VecDoub ptry(ndim);
		double fac1=(1.0-fac)/ndim;
		double fac2=fac1-fac;

		for (int j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

		double ytry=func(ptry);
		if (ytry < y[ihi])
		{
			y[ihi]=ytry;

			for (int j=0;j<ndim;j++)
			{
				psum[j] += ptry[j]-p[ihi][j];
				p[ihi][j]=ptry[j];
			}
		}
		return ytry;
	}
};

#endif
