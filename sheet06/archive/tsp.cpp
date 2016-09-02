#include <iostream>
#include "helfer.h"
#include "anneal.h"
#include "ran.h"

using namespace std;

int main ()
{
	// Anneal and Random object instantiations
	Anneal A;
	Ran R;
	
	// arrays for VecDoub initialization
	double x[] = {3.,2.,0.,3.,1.,0.,2.,1.};
	double y[] = {2.,0.,2.,1.,3.,1.,3.,0.};
	
	// order array
	int r [] = {7,6,5,4,3,2,1,0};

	//int r[] = {0,1,4,5,2,3,6,7};
	
	
	// VecDoub instantiations
	VecDoub xv(8,x);
	VecDoub yv(8,y);
	VecInt rv(8,r);

	// Anneal::order function call
	double dist = A.order(xv,yv,rv);
	
	/* Some remarks:
	 * you can find attached the figure with the positions of the cities and the minimum path found by the program.
	 * However, even though the distance it calculates is correct, the order of the cities is not, since it does
	 * something strange between 7-5 and 2-4. We haven't been able to solve that, but still the numerical distance
	 * for the minimum path length is correct (analitically calculated), that's why we carried on to the next exercise.*/
	
	cout << dist << endl;
	
	cout << "Optimal order:\n";
	for (int i = 0; i < rv.size(); ++i)
		cout << "City: " << rv[i] << endl;
}
