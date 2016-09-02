#include <iostream>
#include "helfer.h"
#include "anneal.h"
#include "ran.h"
#include <ctime>

using namespace std;

int main ()
{
	//Object declarations
	Anneal A;
	Ran R;
	
	// number of cities
	int cities = 1000;
	
	// empty arrays
	double x[cities];
	double y[cities];
	int r[cities];
	
	for (int i = 0; i < cities; ++i)
	{
		// random x and y values
		x[i] = R.doub(0.,10.);
		y[i] = R.doub(0.,5.);
		r[i] = i;
	}
	
	//VecDoub initializations
	VecDoub xv(cities,x);
	VecDoub yv(cities,y);
	VecInt rv(cities,r);
	
	// timers and Anneal::order function call
	clock_t t1 = clock();	
	double dist = A.order(xv,yv,rv);
	clock_t t2 = clock();
	
	// distance and run time
	cout << dist << endl;
	cout << (t2-t1)/double(CLOCKS_PER_SEC) << " seconds"<< endl;
	
	// For 1000 cities:
	// With -O1 compiler optimization: over 30 seconds
	// With -O2: around 17 seconds
	// With -O3: around 14 seconds

}