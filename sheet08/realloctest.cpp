#include <iostream>
#include <cstdlib>	// realloc


using namespace std;


int main () 
{
	int size = 2;
	double* temp = new double[size];
	double* p;
	
	while (temp)
	{
		p = temp;
		temp = (double*)realloc(p,size*size);
		size = size*size;
	}
	
	std::cout << p << std::endl;
	
	
	
}