#include <iostream>
#include "Matrix.h"
#include "Vector.h"
#include "helfer.h"
#include "datum.h"
#include "dynMatrix.h"
#include <ctime>

#include<exception>

using namespace std;

int main() {

	double arr[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
	double arr2[3][3] = {{2,5,6},{1,3,2},{5,3,2}};
	
	Matrix<double> M(arr);
	Matrix<double> D(arr2);
	
	
	
	Vector<double> c(10);
	for (int i=0; i < 10; ++i)
	{
		double t = i;
		c[i] = t;
	}
	
	Vector<double> f = c;
	
	
	int a = 5, b = 2;
	a == MAX(a,b) ? cout << 'a' << endl:cout << 'b' << endl;
	SWAP(a,b);
	
	
	Datum d1;
	cout << d1 << endl;
	
	try 
	{
		cout << d1.weekday() << endl;
	} 
	catch (const char* e)
	{
		cout << e << endl;
	}

}
