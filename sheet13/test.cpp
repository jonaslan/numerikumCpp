// c++11 test program

#include <iostream>
#include <cmath>

using namespace std;

template <typename T, typename S, typename W> 
auto lyapunov(T c, S x, W N) -> decltype(c)
{
	auto sum = 0.;
	for (auto i = 0; i < 10000; ++i)
	{
		sum += log(abs(c*(1.-x[i])));
	}
	return sum/N;
}


int main ()
{
	for (auto i = 0; i < 10000; ++i)
	{
		auto c = 1 + 3./9999.*i;
		auto x = 1e-10;
		for(auto j = 0; j < 10000; ++j)
		{
			x = c*x*(1. - x);
			std::cout << x << std::endl;
		}
	}

	
}
