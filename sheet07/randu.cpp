#include <iostream>
#include <fstream>

int lcg(int a, int c, int x, int m)
{
	return (a*x + c)%m;
}

int main ()
{

	int ai = 4;
	int ci = 0;
	int mi = 200;
	int xi = 1;
	
	for (int i = 0; i < 30; ++i)
	{
		xi = lcg(ai,ci,xi,mi);
		std::cout << xi << std::endl;
	}


	std::ofstream file;
	file.open("randudata.txt",std::ios::trunc);
	int a = 65539;
	int c = 0;
	int m = 2e31;
	int x = 1;
	
	for (int i = 0; i < 10000; ++i)
	{
		x = lcg(a,c,x,m);
		file << x << "\t";
		x = lcg(a,c,x,m);
		file << x << "\t";
		x = lcg(a,c,x,m);
		file << x << "\n";
	}
	file.close();
}