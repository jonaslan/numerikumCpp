#include <iostream>
#include "helfer.h"
#include<cmath>
#include<fstream>
#include "bessel.h"

using namespace std;


int main(){

	std::ofstream file;
	file.open("besseltable.txt",std::ios::trunc);
	
	// two bessel objects
	Bessjy bjy;
	Bessik bik;
	
	double lhs1, rhs1, lhs2, rhs2,lhs3,rhs3, diff1, diff2, diff3;
	file << "#N\tx\tleft side:\tright side:\tdifference:" << std::endl;
	for (int n = 0; n <= 5; ++n)
	{
		double x = 0.;
		
		for (int i = 0; i < 10; ++i)
		{
			
			// calculating the left hand and right hand sides
			
			lhs1 = 2*n*bjy.jn(n,x);
			rhs1 = x*(bjy.jn(n-1,x)+bjy.jn(n+1,x));
			diff1 = lhs1-rhs1;
			lhs2 = 2*n*bik.kn(n,x);
			rhs2 = x*(bik.kn(n+1,x)-bik.kn(n-1,x));
			diff2 = lhs2-rhs2;
			lhs3 = bjy.jn(-n,x);
			rhs3 = pow(-1,n)*bjy.jn(n,x);
			diff3 = lhs3-rhs3;
			
			// Printing the values
			file << "n: " << n << "\t" << "x: " << x << "\n" << "First:\t" << lhs1 << "\t" << rhs1 << "\t" << diff1 << std::endl;
			file << "Second:\t" << lhs2 << "\t" << rhs2 << "\t" << diff2 << std::endl;
			file << "Third:\t" << lhs3 << "\t" << rhs3 << "\t" << diff3 << std::endl << std::endl;
			x += 1.;
		}	
				
		
	}
	file.close();	
	
	
	
	
}