#include <iostream>
#include "Matrix.h"
#include "dynMatrix.h"

using namespace std;

template <typename T>
class Vector
{
	T* v;
	int l;
public:
	Vector(int n)
	{
		v = new T[n];
		l = n;
	}

	Vector(T* ar, int n)
	{
		l = n;
		v = new T[n];
		for (int i = 0; i < n; ++i)
		{
			v[i] = ar[i];
		}
	}

	~Vector()
	{
		delete[] v;
	}

	Vector(const Vector& rhs)
	{
		l = rhs.l;
		v = new T[rhs.l];
	}

	Vector operator+ (Vector& a)
	{
		T* arr = new T[l];
		for (int i = 0; i < l; ++i)
		{
			arr[i] = a.v[i]+v[i];
		}
		Vector retv(arr,l);
		return retv;
	}

	void print() {
		std::cout << "[ ";
		for (int i = 0; i < l; ++i)
		{
			std::cout << v[i] << " ";
		}
		std::cout << "]" << std::endl;
	}
};


int main() {

	double a[3][3] = { {3,0,1},{0,2,4},{1,5,2} };
	double b[3][3] = { {1,0,2},{2,3,0},{1,2,0} };
	Matrix<double> M(a);
	Matrix<double> V(b);

	int arr1[] = {1,2,3};
	int arr2[] = {2,1,4};
	Vector<int> vec1(arr1,3);
	Vector<int> vec2(arr2,3);
	Vector<int> vec3 = vec1+vec2;

	int** pp = new int*[2];
	pp[0] = new int[3];
	pp[1] = new int[3];

	for (int i = 0; i < 2; ++i)
	{
		pp[0][i] = arr1[i];
		pp[1][i] = arr2[i];
	}

	for (int i = 0; i < 2; ++i)
	{
		cout << pp[i][i] << endl;
	}

	int arr3 [] = {1,2,3,4,5,6,7,8,9};
	dynMatrix<int> dM(arr3,3,3);

	//cout << M+V << endl;
	cout << dM+dM;

}
