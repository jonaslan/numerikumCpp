#include <iostream>
#include <cmath>
#include <cstdlib>
#include "minimization.h"
using namespace std;

// Distribution functor
class functor{
	double p,q,v,w;
public:
	functor()
	{
		p = 1.;
		q = 2.;
		v = 100.;
		w = 1.;
	}
	double operator()(const double& x)
	{
		return exp(-(x-p)*(x-p)/v) + exp(-(x+q)*(x+q)/w);
	}
	double df(double& x)
	{
		return((-2./v*(x-p))*exp(-(x-p)*(x-p)/v) + (-2./w*(x+q))*exp(-(x+q)*(x+q)/w));
	}
};




int main ()
{
	functor func;
	Bracketmethod BM;
	
	
	
	// starting values for bracketing
	double x = -0.25;
	double y = 0.25;
	double c;
	Golden G;
	DBrent D;
	Brent B;
	
	BM.bracket(x,y,c,func);
	G.set(x,y,c);
	B.set(x,y,c);
	D.set(x,y,c);
	
	
	double xout = 1;
	double fout;
	int giter = G.minimize(xout,fout,func);
	cout << "Iterations (GOLDEN): " << giter << endl << xout << " " << fout << endl;

	int biter = B.minimize(xout,fout,func);
	cout << "Iterations (BRENT): " << biter << endl << xout << " " << fout << endl;

	int diter = D.minimize(xout, fout,func);
	cout << "Iterations (DBRENT): " << diter << endl << xout << " " << fout << endl;
	
}