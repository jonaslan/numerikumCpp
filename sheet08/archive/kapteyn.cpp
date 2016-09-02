#include <iostream>
#include <cstring>
#include <cmath>
#include "bessel.h"
#include "wynn.h"

using namespace std;

inline double power(const double& x, const int& power)
{
    double r = x;
        for (int i = 1; i < power; ++i)
            r *= x;
        return r;
}

class KapteynSeries
{
private:
    double* values;
    int length;
    int p;
    double x;
    Bessjy BJY;

    void broaden()
    {
        double* temp = new double[length];
        memcpy(temp, values,(sizeof(double)*length));
        delete[] values;
        length = 2*length;
        values = new double[length];
        memcpy(values, temp,(sizeof(double)*length));
        delete[] temp;
    }

    void insert(double input,int n)
    {
        if (p < length)
        {
            values[p] = ::power(n,6)*input + values[p-1];
            ++p;
        } else {
            broaden();
            insert(input,n);
        }
    }

public:
    KapteynSeries(const double& xin, int terms=2) :x(xin),p(0),length(terms)
    {
        Bessjy BJY;
        values = new double[terms];
        this->calculateSum(terms);
    }

    ~KapteynSeries()
    {
        delete[] values;
    }

    double calculateSum(const int& term)
    {
        if (term >= p)
        {
            for (int n = p; n <= term; ++n)
                this->insert( (BJY.jn(n,n*x)*BJY.jn(n,n*x)) ,n);
            return values[term];
        }
        else
        {
            return values[term];
        }
    }
};


// analytic solution:

double analytic(const double& x)
{
    return (x*x*(512. + 27776. * x * x + 161152. * power(x,4) + 189040. * power(x,6) + 45820. * power(x,8) + 1125. * power(x,10)))
    / (2048.*pow((1-x*x),(19./2.)));
}

int main()
{
    const double x =0.5;
    KapteynSeries K(x);


    std::cout.precision(20);

    double diff = 1.;
    const double a = analytic(x);
    int i = 1;

    while (diff > 1e-8)
    {
        diff = std::abs(K.calculateSum(i)-a);
        ++i;
    }


    std::cout << "Sum number " << i << " has value: " << K.calculateSum(i) << std::endl;
    std::cout << "Analytic solution has value: " << analytic(x)<< std::endl;


    /// PART B

    const double epsilon = 10e-8;
    const int nmax = 3000;
    int n = 0;

    Wynn w(nmax,epsilon);

    while (!w.conv && n < nmax)
    {
        ++n;
        double part = K.calculateSum(n);
        double wynn = w.next(part);
    }
    std::cout << "The Wynn-Epsilon Method converges in " <<  n << " iterations" << std::endl;


}
