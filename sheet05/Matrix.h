// if _MATRIX_H_ is not defined, define _MATRIX_H_

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>

// Generalized 3x3 matrix class
template <class T>
class Matrix
{
	// data storage
	T v[3][3];

public:

	// Zero constructor, takes no argument, returns empty matrix
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

	// Constructor, takes multidimensional array of type T as argument
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

	// bracket operator overloading, returns pointer to row i in matrix
	// this operator is used for non-constant access of matrix data, such as assignment
	// and reassignment
	T* operator[] (const int i)
	{
		 return v[i];	
	}
	
	// another bracket operator overload, but this time constant for 
	// constant access to matrix data, such as printing, adding, multiplying
	const T* operator[] (const int i) const
	{
		 return v[i];	
	}
	
	// <<-operator overload for matrix formatting
	friend std::ostream& operator<< (std::ostream& out, const Matrix& mat)
	{
		for (int i = 0; i < 3; ++i)
		{
			out << "[" << mat[i][0] << ",";
			for (int j = 1; j < 2; ++j)
			{
				out << " " << mat[i][j] << ",";
			}
			out << " " << mat[i][2] << "]\n";
		}
		return out;
	}

	// addition operator overload, creates multidimensional array of type T
	// and assigns to it added values of the matrices, returns Matrix
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

	// Multiplication operator overload, multiplies all rows with all columns, and adds
	// values into a multidimensional array. Returns matrix created with same array
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


#endif 