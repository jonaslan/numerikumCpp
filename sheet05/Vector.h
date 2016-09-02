// if _VECTOR_H_ is not defined, define _VECTOR_H_

#ifndef _VECTOR_H_
#define _VECTOR_H_


// Generalized vector class
template <typename T>
class Vector
{
	// Pointer to data array and length of same
	T* data;
	int l;
public:
	// Constructor which only receives size of data array
	Vector(const int N)
	{
		l = N;
		// reservation of heap memory
		data = new T[N];
	}
	
	// Constructor, if undefined, data is set to NULL and length to zero
	Vector(T* ar=NULL, const int N=0)
	{
		// reservation of heap memory
		data = new T[N];
		l = N;
		
		for (int i = 0; i < N; ++i)
		{
			data[i] = ar[i];
		}
	}
	
	// Destructor which deletes data array, thereby freeing allocated memory.
	// The destructor is automatically called whenever a vector object goes out of scope
	// but can also be called manually when data is to be freed
	~Vector()
	{
		delete[] data;
	}

	// copy constructor
	Vector(const Vector& rhs)
	{
		// copy length and reserve new memory
		l = rhs.l;
		data = new T[rhs.l];
		
		// copy content of array
		for (int i = 0; i < rhs.l; ++i)
		{
			data[i] = rhs.data[i];
		}
	}
	
	// bracket operator overloading, returns reference to data 
	T& operator[] (const int& i)
	{
		 return data[i];	
	}

	// addition operator overloading, creates new vector and assigns added values to it
	Vector operator+ (const Vector& a)
	{
		// memory allocation
		T* arr = new T[l];
		
		// add values
		for (int i = 0; i < l; ++i)
		{
			arr[i] = a.data[i]+data[i];
		}
		
		// call to constructor
		Vector R(arr,l);
		return R;
	}

	// <<-operator overloading, optimized for appearence
	friend std::ostream& operator<< (std::ostream& out, Vector v)
	{
		out << "[" << v.data[0] << ",";
		for (int i = 1; i < v.l-1; ++i)
		{
			out << " " << v.data[i] << ",";
		}
		out << " " << v.data[v.l-1] << "]\n";
		return out;
	}
};

#endif