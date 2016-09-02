#include <iostream>
#include <cmath>
#include "helfer.h"
#include "amoeba.h"

using namespace std;


class Functor
{
public:
	VecDoub operator()(VecDoub x)
	{
		return VecDoub(1,-cos(M_PI*(x[0]-1.))-cos(M_PI*(x[1]-1.)));
	}
};

int main ()
{
	Amoeba A;
	Functor func;
	double ar[] = {0.5,0.5};
	VecDoub startV(2,ar);
	double del = 0.1;
	VecDoub point(1,func(startV));
	double fvalmin;
	
	VecDoub minimum = A.minimize(point,del,fvalmin,func);
	
	//cout << point << endl;
}