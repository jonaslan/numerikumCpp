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

	dynMatrix(T* arrin, int rin, int cin)
	{
		v = new T*[r];
		r = rin;
		c = cin;
		int count = 0;
		for (int i = 0; i < r; ++i)
		{
			v[i] = new T[c];
			for (int j = 0; j < c; ++j)
			{
				v[i][j] = arrin[count++];
			}
		}

	}

	~dynMatrix()
	{
		for (int i = 0; i < r; ++i)
		{
			delete[] v[i];
		}
	}

	const T* operator[] (int i) const
	{
		const T* p = v[i];
		return p;
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
		T ret[a.r*a.c];
		int count = 0;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				ret[count++] = v[i][j] + a[i][j];
			}
		}
		dynMatrix R(ret,a.r,a.c);
		return R;
	}

	dynMatrix operator* (const dynMatrix& a)
	{
	}
};

#endif /* DYNMATRIX_H_ */
