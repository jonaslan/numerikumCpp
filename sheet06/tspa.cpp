#include <iostream>
#include "helfer.h"
#include "anneal.h"
#include "ran.h"

// Traveling Salesman Problem, part A

using namespace std;

int main ()
{
	// Anneal and Random object instantiations
	Anneal A;
	Ran R;
	
	// arrays for VecDoub initialization
	double x[] = {1.,2.,0.,1.,3.,0.,2.,3.};
	double y[] = {0.,3.,1.,3.,1.,2.,0.,2.};
	
	// order array
	int r[] = {0,1,2,3,4,5,6,7};
	
	// VecDoub instantiations
	VecDoub xv(8,x);
	VecDoub yv(8,y);
	VecInt rv(8,r);

	// Anneal::order function call
	double path = A.order(xv,yv,rv);
	
	cout << "Path length: " << path << endl;
	
	cout << "Optimal order:\n";
	for (int i = 0; i < rv.size(); ++i)
		cout << "City: " << rv[i] << endl;
	
	
	// The path includes the return to the origin, even though the plot does not.
}
