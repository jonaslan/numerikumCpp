// complex value root finder

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class Z
{

	
	public:
	
	double a;
	double b;
	
	// Constructors for 0, real numbers (0 imaginary part) and complex numbers 
	Z()
	{
		a = 0.;
		b = 0.;
	}
	
	Z(double ain)
	{
		a = ain;
		b = 0.;
	}
	
	Z(double ain, double bin)
	{
		a = ain;
		b = bin;
	}
	
	
	// conjugate function
	Z conjugate(const Z& c)
	{
		return Z(c.a,-c.b);
	}
	
	// addition operator overload
	Z operator +(const Z& in)
	{
		return Z(a+in.a, b+in.b);
	}

	// subtraction operator overload
	Z operator -(const Z& in)
	{
		return Z(a-in.a, b-in.b);
	}
	
	// multiplication operator overload
	Z operator *(const Z& in)
	{
		return Z(a*in.a-b*in.b, a*in.b+b*in.a);
	}
	
	// division operator overload
	Z operator /(const Z& in)
	{
		Z con = conjugate(in);
		Z denom = con*in;
		
		// 'this' points to object in front of operator, '*' dereferences the pointer
		Z numer = *this*con;
			
		return Z(numer.a/denom.a, numer.b/denom.a);
	}
	
	// print function for complex output
	friend ostream& operator<< (ostream& out, const Z c)
	{
		if (c.b >= 0)
		{
			out << c.a << "+" << c.b << "i";
		}
		else
		{
			out << c.a << c.b << "i";
		}
		return out;
	}
};


int main ()
{
	// stepsize for 100 steps from -1 to 1
	double stepsize = 0.02;
	
	// file object
	ofstream data;
	data.open("data.txt",ios::trunc);
	
	// x loop
	for (int i=0; i < 100; ++i)
	{
		double x = -1+i*stepsize;
		
		// y loop
		for (int j=0; j < 100; ++j)
		{
			double y = -1+j*stepsize;
			
			Z z(x,y);
	
			// Newton's method
			for (int k = 0; k < 100; ++k)
			{
				z = z - (z*z*z-1)/(z*z*3);
			}
			
			// check if real part of z is 1
			if (abs(z.a-1.) < 1e-10)
			{
				data << x << "\t" << y << "\t1\n";
			}
			else
			{
				data << x << "\t" << y << "\t0 !!!!!!!!!!!!!!\n";
			}
		}
	}
	data.close();
}
