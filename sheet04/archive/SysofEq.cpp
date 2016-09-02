#include <iostream>

using namespace std;

double* jacobian(double A[])
{
	double* F = new double[4];
	for (int i = 0; i < 4; ++i)
	{
		F[i] = A[i];
	}
	double c = 1./(F[0]*F[3]-F[2]*F[1]);
	F[2] = -F[2];
	F[1] = -F[1];
	for (int i = 0; i < 4; ++i)
	{
		F[i] = F[i]*c;
	}
	return F;
}

int main ()
{
	double P[] = {1,2,3,2};		// really 2X2
	double* J = jacobian(P);
	cout << J[1];
	
	delete J;
}
