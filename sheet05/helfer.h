#ifndef _HELFER_H_
#define _HELFER_H_

// Generalized max,sign, swap functions
template <typename T>
T MAX(const T& a,const T& b)
{
	// ternary operator ? evaluates a>b, if true it returns a, if false b
	return a>b ? a:b;
}

template <typename T>

// SWAP function uses a temporary variable of type T for storing a, then replaces a with b and
// b with temp
void SWAP(T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;
}


template <typename T>

// SIGN function, evaluates four conditions and returns a with appropriate sign
T SIGN(const T& a, const T& b)
{
	return a >= 0 && b >= 0 ? a : a < 0 && b >= 0 ? -a : a >= 0 && b < 0 ? -a : a;
}

#endif