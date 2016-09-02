#ifndef _FFT_H_
#define _FFT_H_

#include <cmath>
#include <complex>
#include "helfer.h"
using namespace std;

typedef complex<double> doz;
typedef Vector<doz> VecDoz;

// **********
void fft(double data[], const int n, const int isign)
{
	if (n<2 || n&(n-1))
		throw("Number of data points must be a power of 2!");

	int nn = n << 1;
	int j = 1, s;

	for (int i=1; i<nn; i+=2)
	{
		if (j > i)
		{
			SWAP(data[j-1], data[i-1]);
			SWAP(data[j], data[i]);
		}

		s = n;
		while (s >= 2 && j > s)
		{
			j -= s;
			s >>= 1;
		}
		j += s;
	}

	int istep, mmax = 2;
	double delta, wtemp, tempr, tempi, alpha, beta, wr, wi;

	while (mmax < nn)
	{
		istep = mmax << 1;
		delta = 2. * M_PI * isign / mmax;

		wtemp = sin(0.5 * delta);
		alpha = 2. * wtemp*wtemp;
		beta = sin(delta);

		wr = 1.;
		wi = 0.;

		for (int m=1; m<mmax; m+=2)
		{
			for (int i=m; i<=nn; i+=istep)
			{
				int k = i + mmax;
				tempr = wr * data[k-1] - wi * data[k];
				tempi = wr * data[k]   + wi * data[k-1];

				data[k-1] = data[i-1] - tempr;
				data[k]   = data[i]   - tempi;

				data[i-1] += tempr;
				data[i]   += tempi;
			}

			wtemp = wr;
			wr = wr - ( alpha * wtemp + beta * wi    );
			wi = wi - ( alpha * wi    - beta * wtemp );
		}

		mmax = istep;
	}
	
}

// **********
void fft(VecDoub& data, const int isign)
{
	fft(&data[0], data.size()/2, isign);
}

// **********
void fft(VecDoz& data, const int isign)
{
	const int n = data.size();
	double *tmp = new double[2*n];

	for (int i=0; i<n; ++i)
	{
		tmp[2*i]   = data[i].real();
		tmp[2*i+1] = data[i].imag();
	}

	fft(tmp, n, isign);

	for (int i=0; i<n; ++i)
	{
		double real = tmp[2*i];
		double imag = tmp[2*i+1];
		data[i] = doz(real,imag);
	}
	delete[] tmp;
}

// **********
void realfft(VecDoub& data, const int isign)
{
	int i, i1, i2, i3, i4;
	int n = data.size();

	double c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp;
	double theta = M_PI/double(n>>1);

	if (isign == 1)
	{
		c2 = -0.5;
		fft(data,1);
	}

	else
	{
		c2 = 0.5;
		theta = -theta;
	}

	wtemp = sin(0.5*theta);
	wpr = -2.0 *wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;

	for (i=1; i<(n>>2); i++)
	{
		i2 = 1 + (i1=i+i);
		i4 = 1 + (i3=n-i1);

		h1r = c1 * ( data[i1]+data[i3] );
		h1i = c1 * ( data[i2]-data[i4] );
		h2r = -c2 * ( data[i2]+data[i4] );
		h2i = c2 * ( data[i1]-data[i3] );

		data[i1] = h1r + wr*h2r - wi*h2i;
		data[i2] = h1i + wr*h2i + wi*h2r;
		data[i3] = h1r - wr*h2r + wi*h2i;
		data[i4] = -h1i + wr*h2i + wi*h2r;

		wr = (wtemp=wr)*wpr - wi*wpi + wr;
		wi = wi*wpr + wtemp*wpi + wi;
	}

	if (isign == 1)
	{
		data[0] = (h1r=data[0]) + data[1];
		data[1] = h1r - data[1];
	}

	else
	{
		data[0] = c1 * ( (h1r=data[0])+data[1] );
		data[1] = c1 * (h1r-data[1]);

		fft(data,-1);
	}
}

#endif
