#ifndef _FEM_1D_
#define _FEM_1D_

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "helfer.h"
#include "adapt.h"
#include "ludcmp.h"
#include "integrand.h"

enum { SLAB, ISO };

// **********
template <class T>
class Fem1D
{
	T &f;
	const int Nx, Mt;
	const double dt;
	const double acc;
	
	VecDoub coordinates, u;
	MatInt elements;

	double dx;
	int digits;
	ofstream *datei;

// 	Initialization and file output
	void uInit(void);
	void fileOutput(const int tIdx);

// 	Calculate next time step
	void step(const int tIdx);

// 	Coefficients and integrands
	double coeffA(const VecDoub& x, const int j, const int k);
	double coeffB(const VecDoub& x, const int j, const int k);
	double coeffC(const VecDoub& x, const int j, const int k);
	double coeffF(const VecDoub& x, const int k);

	public:
// 	Constructor: functor, number of spatial/time steps, time increment dt, accuracy
	Fem1D(T& func, const int NxIn, const int MtIn, const double dtIn, const double accIn=1e-10);
	~Fem1D(void);

// 	Run integration: return vectors for norm and mean-square displacement
// 	Definition of the norm(t): integral from -1 to 1 of f(t,mu) over mu
// 	Definition of the msq(t): integral from -1 to 1 of (mu-mu0)^2*f(t,mu) over mu
	void run(VecDoub& nrm, VecDoub& msq);
};

// **********
template <class T>
Fem1D<T>::Fem1D(T& func, const int NxIn, const int MtIn, const double dtIn, const double accIn) : f(func), \
	Nx(NxIn), Mt(MtIn), dt(dtIn), acc(accIn), coordinates(Nx), u(Nx), elements(Nx-1,2)
{
// 	Spatial domain and increment
	const double a = -1.;
	const double b = 1.;
	dx = (b-a)/(Nx-1);

//	Coordinates
	for (int i=0; i<Nx; ++i)
		coordinates[i] = a + i*dx;

//	Elements
	for (int i=0; i<Nx-1; ++i)
		for (int j=0; j<2; ++j)
			elements[i][j] = i+j;

// 	How many digits should the file names have?
	digits = static_cast<int>(log10(Mt)+1);
	datei = new ofstream[Mt];

// 	Initial function u(x,t=0)
	uInit();
}

// **********
template <class T>
Fem1D<T>::~Fem1D()
{
	delete[] datei;
}

// **********
template <class T>
void Fem1D<T>::run(VecDoub& nrm, VecDoub& msq)
{
	nrm.assign(Mt,0.);
	msq.assign(Mt,0.);

	for (int tIdx=0; tIdx<Mt; ++tIdx)
	{
		// Calculate next time step
		step(tIdx);

		// Calculate norm and mean-square displacement
		for (int j=0; j<Nx-1; ++j)
		{
			nrm[tIdx] += dx/2. * ( u[elements[j][0]]+u[elements[j][1]] );
			msq[tIdx] += dx/2. * ( u[elements[j][0]]+u[elements[j][1]] ) * SQR( coordinates[elements[j][0]]+coordinates[elements[j][1]] - 1. )/4.;
		}

		// Print time step to file
		fileOutput(tIdx);
	}
}

// **********
template <class T>
void Fem1D<T>::step(const int tIdx)
{
	MatDoub A(Nx,Nx,0.), B(Nx,Nx,0.), C(Nx,Nx,0.);
	VecDoub b(Nx,0.);

//	Assembly
	for (int i=0; i<Nx-1; ++i)
	{
		VecDoub x(2);
		for (int j=0; j<2; ++j)
			x[j] = coordinates[elements[i][j]];
		
		for (int j=0; j<2; ++j)
			for (int k=0; k<2; ++k)
			{
				A[elements[i][j]][elements[i][k]] += coeffA(x,j,k);
				B[elements[i][j]][elements[i][k]] += coeffB(x,j,k);
				C[elements[i][j]][elements[i][k]] += coeffC(x,j,k);
			}

		for (int k=0; k<2; ++k)
			b[elements[i][k]] += coeffF(x,k);
	}
	
// 	Calculate b = C * u_old
	for (int i=0; i<Nx; ++i)
	{
		double sum = 0.;
		for (int j=0; j<Nx; ++j)
			sum += C[i][j] * u[j];
		b[i] += sum;
	}

// 	Calculate M = A*dt + C
	MatDoub M(Nx,Nx,0.);
	for (int i=0; i<Nx; ++i)
		for (int j=0; j<Nx; ++j)
			M[i][j] = A[i][j]*dt + C[i][j];

//	Dirichlet conditions
	u.assign(Nx,0.);

// 	Calculate b = b - (A*dt+C) * u;
	for (int i=0; i<Nx; ++i)
	{
		double sum = 0.;
		for (int j=0; j<Nx; ++j)
			sum += (A[i][j]*dt + C[i][j]) * u[j];
		b[i] -= sum;
	}

	VecDoub g(Nx,1.);
	g[0] = g[Nx-1] = 0.5;

// 	Generate extended matrix MM = ([M;g]^T * [M;g])
	MatDoub MM(Nx,Nx,0.);
	for (int i=0; i<Nx; ++i)
		for (int j=0; j<Nx; ++j)
			for (int k=0; k<Nx+1; ++k)
				MM[i][j] += k==Nx ? g[i]*g[j]*SQR(dx) : M[k][i]*M[k][j];

// 	Generate extended RHS vector bb = ([M;g]^T * b)
	VecDoub bb(Nx,0.);
	for (int i=0; i<Nx; ++i)
		for (int k=0; k<Nx+1; ++k)
			bb[i] += k==Nx ? g[i]*dx : M[k][i]*b[k];

// 	Solve overdetermined system MM*u = bb
	LUdcmp lu(MM);
	lu.solve(bb,u);
	lu.mprove(bb,u);

	double norm = 0.;
	for (int i=0; i<Nx; ++i)
		norm += g[i]*dx * u[i];

	std::cout << "Finished step #" << tIdx+1 << " of " << Mt << " total steps.";
	std::cout << std::endl;
}

// **********
template <class T>
double Fem1D<T>::coeffA(const VecDoub& x, const int j, const int k)
{
	const int vz1 = pow(-1,j);
	const int vz2 = pow(-1,k);

// 	Spawn integrand with coefficient function a
	Adapt adapt(acc);
	Integrand0<T> i( &f, &T::a );

	double integral = adapt.integrate( i, x[0], x[1] );
	integral /= SQR(x[1]-x[0]);

	return vz1 * vz2 * integral;
}

// **********
template <class T>
double Fem1D<T>::coeffB(const VecDoub& x, const int j, const int k)
{
	const int vz1 = pow(-1,j+1);
	const int vz2 = pow(-1,k+1);

// 	Spawn integrand with coefficient function b
	Adapt adapt(acc);
	Integrand1<T> i( &f, &T::b, 1., -x[1-j%2], 0., 1. );

	double integral = adapt.integrate( i, x[0], x[1] );
	integral /= SQR(x[1]-x[0]);

	return vz1 * vz2 * integral;
}

// **********
template <class T>
double Fem1D<T>::coeffC(const VecDoub& x, const int j, const int k)
{
	const int vz1 = pow(-1,j+1);
	const int vz2 = pow(-1,k+1);

// 	Spawn integrand with coefficient function c
	Adapt adapt(acc);
	Integrand1<T> i( &f, &T::c, 1., -x[1-j%2], 1., -x[1-k%2] );

	double integral = adapt.integrate( i, x[0], x[1] );
	integral /= SQR(x[1]-x[0]);

	return vz1 * vz2 * integral;
}

// **********
template <class T>
double Fem1D<T>::coeffF(const VecDoub& x, const int k)
{
	int vz = pow(-1,k+1);

// 	Spawn integrand with coefficient function f
	Adapt adapt(acc);
	Integrand1<T> i( &f, &T::f, 0., 1., 1., -x[1-k%2] );

	double integral = adapt.integrate(i,x[0],x[1]);
	integral /= x[1]-x[0];

	return vz * integral;
}

// **********
template <class T>
void Fem1D<T>::uInit()
{
//	Volume force in considered domain
	int uIdx = Nx/2 + f.m0/dx;
	if (uIdx == Nx) --uIdx;
	else if (uIdx < 0) ++uIdx;

	u[uIdx] = 1./dx;
}

// **********
template <class T>
void Fem1D<T>::fileOutput(const int tIdx)
{
	stringstream filen;
	filen << "step";
	filen << std::setfill('0') << setw(digits) << tIdx+1;
	filen << ".txt";

	datei[tIdx].open(filen.str().c_str(),std::ios::out|std::ios::trunc);
	for (int j=0; j<Nx-1; ++j)
		datei[tIdx] << coordinates[elements[j][0]] << '\t' << u[elements[j][0]] << std::endl;
	datei[tIdx] << coordinates[elements[Nx-2][1]] << '\t' << u[elements[Nx-2][1]] << std::endl;
	datei[tIdx].close();
}

#endif
