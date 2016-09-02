#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>

template <class T>
class Matrix
{
	T v[3][3];

public:

	Matrix()
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				v[i][j] = 0;
			}
		}
	}

	Matrix(T a[3][3])
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				v[i][j] = a[i][j];
			}
		}
	}

	const T* operator[] (int i) const
	{
		const T* p = v[i];
		return p;
	}

	friend std::ostream& operator<< (std::ostream& out, Matrix mat)
	{
		for (int i = 0; i < 3; ++i)
		{
			out << "[";
			for (int j = 0; j < 3; ++j)
			{
				out << " " << mat[i][j] << ",";
			}
			out << "]\n";
		}
		return out;
	}

	Matrix operator+ (const Matrix& a)
	{
		T ret[3][3];
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				ret[i][j] = v[i][j] + a[i][j];
			}
		}
		Matrix R(ret);
		return R;
	}

	Matrix operator* (const Matrix& a)
	{
		T ret[3][3] = {{0, 0, 0}, {0, 0, 0},{0, 0, 0}};
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 3; ++k)
				{
						ret[i][j] += v[i][k]*a[k][j];
				}
			}
		}
		Matrix R(ret);
		return R;
	}
};





#endif /* MATRIX_H_ */
