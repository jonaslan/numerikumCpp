#include <iostream>
#include <cmath>
#include <cstdlib>
#include "minimization.h"
using namespace std;

// Distribution functor class
class functor{
	double p,q,v,w;
public:
	functor()
	{
		//Initialisation of parameters in f(x)
		p = 1.;
		q = 2.;
		v = 100.;
		w = 1.;
	}
	
	// ()-operator overload
	double operator()(const double& x)
	{
		return exp(-(x-p)*(x-p)/v) + exp(-(x+q)*(x+q)/w);
	}
	
	//Definition of the derivative of f(x), required by Brent's method with the use of derivative
	double df(double& x)
	{
		return((-2./v*(x-p))*exp(-(x-p)*(x-p)/v) + (-2./w*(x+q))*exp(-(x+q)*(x+q)/w));
	}
};




int main ()
{
	//Instantiation of functor
	functor func;
	
	// creation of objects specific to each numerical method
	Golden G;
	DBrent D;
	Brent B;
	Bracketmethod BM;
	
	// starting values for bracketing
	double x = -0.25;
	double y = 0.25;
	double c;
	
	//bracketing
	BM.bracket(x,y,c,func);
	
	//set the initial values provided by bracketing to the different methods
	G.set(x,y,c);
	B.set(x,y,c);
	D.set(x,y,c);
	
	//Initialise variable that will store position of minimum
	double xout = 1;
	
	//Initialise variable that will store value of the minimum
	double fout;
	
	//Call to different methods
	int giter = G.minimize(xout,fout,func);
	int biter = B.minimize(xout,fout,func);
	int diter = D.minimize(xout, fout,func);
	
	//Display the number of iterations required by each method to converge.
	/*We observed the two different versions of the Brent method speed up convergence
	 * the one that used the derivative even more since it has additional information
	 * about the function*/
	cout << "Iterations (GOLDEN): " << giter << endl << xout << " " << fout << endl;
	cout << "Iterations (BRENT): " << biter << endl << xout << " " << fout << endl;
	cout << "Iterations (DBRENT): " << diter << endl << xout << " " << fout << endl;
	
}