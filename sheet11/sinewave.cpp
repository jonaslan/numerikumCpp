#include <iostream>
#include <fstream>
#include "fft.h"
#include "helfer.h"


class Functor
{
private:
	double a1,a2,a3;

public:
	Functor() : a1(0.05),a2(0.2),a3(0.6) { }
	
	double operator () (double x)
	{
		return sin(a1*x) + sin(a2*x) + sin(a3*x);
	}
};

int main () {
	const double xmin = 0;
	const double xmax = 1000;
	Functor f;
	std::ofstream file;
	file.open("sinewave.txt", std::ios::trunc);
	
	
	int n = 9;
	int N = pow(2,n);
	
	double lambda = (xmax-xmin)/N;
	
	VecDoub data(2*N,0.);
	
	VecDoub freq(2*N,0.);
	
	
	
	for (int i = 0; i < 2*N; i+=2)
	{
		data[i] = f(xmin+lambda*i);
		if (i >= 2 && i <= N)
		{
			freq[i] = freq[i-1]+1./(N*lambda);
			freq[i+1] = freq[i];
		}
		
		else if (i > N)
		{
			freq[i] = -freq[2*N-i];
			//freq[i]=-1./(2*lambda)+1/(N*lambda)*
			freq[i+1] = freq[i];
		}
	}
	
	fft(data,1);
	/*
	for (int i = 0; i < data.size(); i+=2)
		std::cout << freq[i] << "\t" << data[i] <<   std::endl;
	
	*/
	//for (int i = 1; i < data.size(); i+=2)
	//	std::cout << freq[i] << "\t" << data[i] << std::endl;
	
	
	/* Print to file
	 *
	 *  freq = 0;
	 *for (int i = 0; i < 2*N; i +=2)
	 {
		 file << data[i] << "\t";
		
		if (i >= 2 && i <= N)
		{
			freq = freq+1./(N*lambda);
		}
		if(i == N)
			freq = -freq;
		else if (i > N)
		{
			freq = freq + 1./(N*lambda);
		}
		
		file << freq << std::endl;
	 *
	 * 
	 */
	
	
	fft(data,-1);
	
	for (int i = 0; i < data.size(); i+=2)
		std::cout << freq[i] << "\t" << data[i] << std::endl;

	
	file.close();
}