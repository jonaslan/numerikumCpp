#ifndef DYNMATRIX_H_
#define DYNMATRIX_H_

#include <iostream>
#include <algorithm>
template <class T>
class dynMatrix
{
	T** v;
	int r;
	int c;

public:


	dynMatrix()
	{
		*v = NULL;
		r = 0;
		c = 0;
	}

	dynMatrix(int rin, int cin)
	{
		v = new T*[r];
		r = rin;
		c = cin;
		//int count = 0;
		for (int i = 0; i < r; ++i)
		{
			v[i] = new T[c];
		}

	}

	~dynMatrix()
	{
		for (int i = 0; i < r; ++i)
		{
			delete[] v[i];
		}
	}

	 T* operator[] (const int i)
	{
		 return v[i];	
	}
	
	const T* operator[] (const int i) const
	{
		 return v[i];	
	}
	
	friend std::ostream& operator<< (std::ostream& out, const dynMatrix& mat)
	{
		for (int i = 0; i < mat.r; ++i)
		{
			out << "[";
			for (int j = 0; j < mat.c; ++j)
			{
				out << " " << mat[i][j] << ",";
			}
			out << "]\n";
		}
		return out;
	}

	dynMatrix operator+ (const dynMatrix& a)
	{
		dynMatrix R(a.r,a.c);
		
		for (int i = 0; i < a.r; ++i)
		{
			for (int j = 0; j < a.c; ++j)
			{
				R[i][j] = v[i][j] + a[i][j];
			}
		}
		
		return R;
	}

	dynMatrix operator* (const dynMatrix& a)
	{
		dynMatrix R(c,a.r);
		if (c == a.r && v != NULL && a.v != NULL)
		{

			for (int i = 0; i < c; ++i)
			{
				for (int j = 0; j < a.r; ++j)
				{
					R[i][j] = 0;
				}
			}
			for (int i = 0; i < c; ++i)
			{
				for (int j = 0; j < a.r; ++j)
				{
					for (int k = 0; k < c; ++k)
					{
						R[i][j] += v[i][k]*a[k][j];
					}
				}
			}
			std::cout << R << std::endl;
			return R;
		}
		else
		{
			throw("Invalid matrix operation");
		}

	}
};

#endif /* DYNMATRIX_H_ */
