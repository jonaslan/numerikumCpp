#include <iostream>
#include <cmath>
#include "helfer.h"
#include "amoeba.h"

using namespace std;


// functor class
class Functor
{
public:
	// ()-operator overload
	double operator()(const VecDoub& x)
	{
		return -cos(M_PI*(x[0]-1.))-cos(M_PI*(x[1]-1.));
	}
};

int main ()
{
	// auxiliary variables
	double del = 0.1;
	double fvalmin;
	
	// Accuracy and starting value inputs
	cout << "Enter accuracy: \n";
	double a;
	cin >> a;
	double ar[2] = {0};
	cout << "Enter starting values:\n";
	cin >> ar[0];
	cin >> ar[1];

	// Amoeba object instantiation with tolerance as passed parameter	
	Amoeba A(a);
	
	// functor instantiation
	Functor func;
	
	// VecDoub of starting values
	VecDoub startV(2,ar);
	
	// Amoeba::minimize function call
	VecDoub minimum = A.minimize(startV,del,fvalmin,func);
	
	cout << minimum[0] << " " << minimum[1] << " " << fvalmin << endl;
}