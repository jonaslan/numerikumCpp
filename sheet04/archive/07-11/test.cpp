// Einklammern und Bisektionsverfahren
 
// N = 100, c = 1,1
 
 
#include <iostream>
#include <cmath> //exp()
#include<cstdio> // printf()
 
using namespace std;
 
static int N = 100;
static double redFactor = 1.1;
 
class Functor
{
    double a;
 
    public:
     
    Functor(double ain)
    {
        a = ain;
    }
    double operator () (double x)
    {
        return exp(x)-a*x;
    }
};
 
 
int bracket (Functor& func, double& a, double& b)
{
    double f = func(a);
    double g = func(b);
     
    for (int i = 0; i < N; ++i)
    {
        if (f*g < 0)
        {
            return i;
        }
        if (abs(f) < abs(g))
        {
            f = func( a+= redFactor*(a-b));
        }
        else
        {
            g = func( b += redFactor*(b-a));
        }
    }
    return 0;
}
 
int root(Functor& myfunc, double& a, double&b, double& ans)
{
    double mid = (a+b)/2.;
    double diff = mid-a;
    int count = 0;
     
    while (diff > 1e-14 && count < 1000)
    {
        if (bracket(myfunc,a,mid))
        {
            b = mid;
            mid = (a+b)/2.;
            diff = mid-a;
        }
        else
        {
            a = mid;
            mid = (a+b)/2.;
            diff = mid-a;
        }
        ++count;
    }
    ans = mid;
    return count;
}
 
int main () {
    Functor myfunc(4.);
    double a = 1.0;
    double b = 1.5;
    double zero;
     
    int iter = bracket(myfunc, a, b);
     
    if (iter == 0)
    {
        cout << "Failed to bracket\n";
    }
    else
    {
        cout << iter << " iterations.\n";
        cout << "a: " << a << " b: " << b << endl;
     
        // bisection
     
    iter = root(myfunc, a, b, zero);
    double fval = myfunc(zero);
 
    printf("Root: %f, Function value at root: %f, reached in: %d iterations", zero, fval, iter);
    }
}
