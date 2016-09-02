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
	
	Z conjugate(const Z& c)
	{
		return Z(c.a,-c.b);
	}
	
	
	Z operator +(const Z& in)
	{
		return Z(a+in.a, b+in.b);
	}

	Z operator -(const Z& in)
	{
		return Z(a-in.a, b-in.b);
	}
	
	Z operator *(const Z& in)
	{
		return Z(a*in.a-b*in.b, a*in.b+b*in.a);
	}
	
	Z operator /(const Z& in)
	{
		Z con = conjugate(in);
		Z denom = con*in;
		Z numer = *this*con;
			
		return Z(numer.a/denom.a, numer.b/denom.a);
	}
	
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
	
	double absVal()
	{
		return sqrt(a*a+b*b);
	}
};


int main ()
{
	double stepsize = 0.02;
	ofstream data;
	data.open("data.txt",ios::trunc);
	
	for (int i=0; i < 100; ++i)
	{
		double x = -1+i*0.02;
		for (int j=0; j < 100; ++j)
		{
			double y = -1+j*0.02;
			
			Z z(x,y);
	
			for (int k = 0; k < 100; ++k)
			{
				z = z - (z*z*z-1)/(z*z*3);
			}
			if ((z.a-1.) < 1e-10)
			{
				data << x << "\t" << y << "\t1\n";
			}
			else
			{
				data << x << "\t" << y << "\t0 !!!!!!!!!!!!!!\n";
			}
		}
	}
}
