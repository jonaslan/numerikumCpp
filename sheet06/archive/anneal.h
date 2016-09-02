#ifndef _ANNEAL_H_
#define _ANNEAL_H_

#include "ran.h"
#include "helfer.h"

// **********
class Anneal
{
	Ran myran;

// 	List of auxiliary functions that are defined outside of the class
	double revcst(VecDoub &x, VecDoub &y, VecInt &iorder, VecInt &n);
	void reverse(VecInt &iorder, VecInt &n);
	double trncst(VecDoub &x, VecDoub &y, VecInt &iorder, VecInt &n);
	void trnspt(VecInt &iorder, VecInt &n);
	bool metrop(const double de, const double t);
	inline double alen(const double a, const double b, const double c, const double d);

	public:
//	Optimization method to be called by user:
//	x and y are VecDoub vectors containing the coordinates;
//	iorder is a VecInt vector in which the optimal sequence will be stored
	double order(VecDoub &x, VecDoub &y, VecInt &iorder);
};

// **********
double Anneal::order(VecDoub &x, VecDoub &y, VecInt &iorder)
{
	// Temperatur zu Beginn
	double t = 0.5;

	// Verringerung der Temperatur pro Schritt
	const double abk = 0.9;

	bool ans;
	int i,i1,i2,nn;
	VecInt n(6);
	double de;
	double path=0.0;

	int ncity = x.size();
	int nover = 1000*ncity;
	int nlimit = 100*ncity;

	for (i=0;i<ncity-1;i++)
	{
		i1 = iorder[i];
		i2 = iorder[i+1];
		path += alen(x[i1],x[i2],y[i1],y[i2]);
	}

	i1 = iorder[ncity-1];
	i2 =  iorder[0];
	path += alen(x[i1],x[i2],y[i1],y[i2]);

	for (int j=0;j<100;j++)
	{
		int nsucc = 0;
		for (int k=0;k<nover;k++)
		{
			do
			{
				n[0] = int(ncity*myran.doub());
				n[1] = int((ncity-1)*myran.doub());

				if (n[1] >= n[0])
					++n[1];

				nn = (n[0]-n[1]+ncity-1) % ncity;
			}
			while (nn<2);

			if (myran.doub() < 0.5)
			{
				n[2] = n[1]+int(abs(nn-1)*myran.doub())+1;
				n[2] %= ncity;

				de = trncst(x,y,iorder,n);
				ans = metrop(de,t);

				if (ans)
				{
					++nsucc;
					path += de;
					trnspt(iorder,n);
				}
			}

			else
			{
				de = revcst(x,y,iorder,n);
				ans = metrop(de,t);

				if (ans)
				{
					++nsucc;
					path += de;
					reverse(iorder,n);
				}
			}

			if (nsucc >= nlimit)
				break;
		}

		t *= abk;
		if (nsucc == 0)
			return path;
	}
}

// **********
double Anneal::revcst(VecDoub &x, VecDoub &y, VecInt &iorder, VecInt &n)
{
	VecDoub xx(4),yy(4);
	int ncity=x.size();
	n[2]=(n[0]+ncity-1) % ncity;
	n[3]=(n[1]+1) % ncity;

	for (int j=0;j<4;j++)
	{
		int ii = iorder[n[j]];
		xx[j] = x[ii];
		yy[j] = y[ii];
	}

	double de = -alen(xx[0],xx[2],yy[0],yy[2]);
	de -= alen(xx[1],xx[3],yy[1],yy[3]);
	de += alen(xx[0],xx[3],yy[0],yy[3]);
	de += alen(xx[1],xx[2],yy[1],yy[2]);
	return de;
}

// **********
void Anneal::reverse(VecInt &iorder, VecInt &n)
{
	int ncity = iorder.size();
	int nn = (1+((n[1]-n[0]+ncity) % ncity))/2;

	for (int j=0;j<nn;j++)
	{
		int k = (n[0]+j) % ncity;
		int l = (n[1]-j+ncity) % ncity;
		int itmp = iorder[k];
		iorder[k] = iorder[l];
		iorder[l] = itmp;
	}
}
	
// **********
double Anneal::trncst(VecDoub &x, VecDoub &y, VecInt &iorder, VecInt &n)
{
	VecDoub xx(6),yy(6);
	int ncity = x.size();

	n[3] = (n[2]+1) % ncity;
	n[4] = (n[0]+ncity-1) % ncity;
	n[5] = (n[1]+1) % ncity;

	for (int j=0;j<6;j++)
	{
		int ii = iorder[n[j]];
		xx[j] = x[ii];
		yy[j] = y[ii];
	}

	double de = -alen(xx[1],xx[5],yy[1],yy[5]);
	de -= alen(xx[0],xx[4],yy[0],yy[4]);
	de -= alen(xx[2],xx[3],yy[2],yy[3]);
	de += alen(xx[0],xx[2],yy[0],yy[2]);
	de += alen(xx[1],xx[3],yy[1],yy[3]);
	de += alen(xx[4],xx[5],yy[4],yy[5]);
	return de;
}

// **********
void Anneal::trnspt(VecInt &iorder, VecInt &n)
{
	int ncity = iorder.size();
	VecInt jorder(ncity);
	int m1 = (n[1]-n[0]+ncity) % ncity;
	int m2 = (n[4]-n[3]+ncity) % ncity;
	int m3 = (n[2]-n[5]+ncity) % ncity;
	int nn=0;

	for (int j=0;j<=m1;j++)
	{
		int jj = (j+n[0]) % ncity;
		jorder[nn++] = iorder[jj];
	}

	for (int j=0;j<=m2;j++)
	{
		int jj = (j+n[3]) % ncity;
		jorder[nn++] = iorder[jj];
	}

	for (int j=0;j<=m3;j++)
	{
		int jj = (j+n[5]) % ncity;
		jorder[nn++] = iorder[jj];
	}

	for (int j=0;j<ncity;j++)
		iorder[j] = jorder[j];
}

// **********
bool Anneal::metrop(const double de, const double t)
{
	return de < 0.0 || myran.doub() < exp(-de/t);
}

// **********
inline double Anneal::alen(const double a, const double b, const double c, const double d)
{
	return sqrt((b-a)*(b-a)+(d-c)*(d-c));
}

#endif
