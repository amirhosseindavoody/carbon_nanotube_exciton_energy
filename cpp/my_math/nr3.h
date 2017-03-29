#ifndef _NR3_H_
#define _NR3_H_

//#define _CHECKBOUNDS_ 1
//#define _USESTDVECTOR_ 1
//#define _USENRERRORCLASS_ 1
//#define _TURNONFPES_ 1

// all the system #include's we'll ever need
#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

// using namespace std;

namespace nr
{

// macro-like inline functions

template<class T>
inline T sqr(const T a) {return a*a;}

template<class T>
inline const T &max(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float max(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float max(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &min(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float min(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float min(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline T sign(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float sign(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float sign(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void swap(T &a, T &b)
	{T dum=a; a=b; b=dum;}

// exception handling

#ifndef _USENRERRORCLASS_
#define throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
#else
struct NRerror {
	char *message;
	char *file;
	int line;
	NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};
#define throw(message) throw(NRerror(message,__FILE__,__LINE__));
void NRcatch(NRerror err) {
	printf("ERROR: %s\n     in file %s at line %d\n",
		err.message, err.file, err.line);
	exit(1);
}
#endif

// usage example:
//
//	try {
//		somebadroutine();
//	}
//	catch(NRerror s) {NRcatch(s);}
//
// (You can of course substitute any other catch body for NRcatch(s).)


// Vector and Matrix Classes

#ifdef _USESTDVECTOR_
#define NRvector std::vector
#else

template <class T>
class NRvector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRvector();
	explicit NRvector(int n);		// Zero-based array
	NRvector(int n, const T &a);	// Initialize to constant value
	NRvector(int n, const T *a);	// Initialize to array
	NRvector(const NRvector &rhs);	// Copy constructor
	NRvector & operator=(const NRvector &rhs);	//Copy assignment
	NRvector(NRvector &&rhs);	// Move constructor
	NRvector & operator=(NRvector &&rhs);	//Move assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	NRvector operator+(const NRvector &other);	// Addition operator
	NRvector operator*(const T &a) const;	// Multiply by number operator
	~NRvector();
};

// NRvector definitions

template <class T>
NRvector<T>::NRvector() : nn(0), v(NULL) {}

// Zero-based array
template <class T>
NRvector<T>::NRvector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

// Initialize to constant value
template <class T>
NRvector<T>::NRvector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = a;
}

// Initialize to array
template <class T>
NRvector<T>::NRvector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}

// Copy constructor
template <class T>
NRvector<T>::NRvector(const NRvector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

// Move constructor
template <class T>
NRvector<T>::NRvector(NRvector<T> &&rhs) : nn(rhs.nn), v(nn>0 ? rhs.v : NULL)
{
	rhs.nn = 0;
	rhs.v = NULL;
}

// Copy assignment
template <class T>
NRvector<T> & NRvector<T>::operator=(const NRvector<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != NULL) delete [] (v);
			nn=rhs.nn;
			v= nn>0 ? new T[nn] : NULL;
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

// Move assignment
template <class T>
NRvector<T> & NRvector<T>::operator=(NRvector<T> &&rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (v != NULL) delete [] (v);
		nn=rhs.nn;
		v= nn>0 ? rhs.v : NULL;
		rhs.v = NULL;
	}
	return *this;
}

//i'th element
template <class T>
inline T & NRvector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

//i'th element
template <class T>
inline const T & NRvector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRvector subscript out of bounds");
}
#endif
	return v[i];
}

// size of the vector
template <class T>
inline int NRvector<T>::size() const
{
	return nn;
}

// resize the vector do not vector values
template <class T>
void NRvector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

// resize and assign a constant value
template <class T>
void NRvector<T>::assign(int newn, const T& a)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
	for (int i=0;i<nn;i++) v[i] = a;
}

// Addition operator
template <class T>
NRvector<T> NRvector<T>::operator+(const NRvector<T> &other)
{
	if (v==NULL) throw("Trying to add empty vector!");
	if (other.v==NULL) throw("Trying to add empty vector!");
	if (nn != other.nn) throw("NRvector dimensions do not match for addition!");

	NRvector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i] + other[i];

	return res;
}

// Multiply by vector by number
template <class T>
NRvector<T> NRvector<T>::operator*(const T &a) const
{
	if (v==NULL) throw("Trying to multiply empty vector!");

	NRvector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = a * v[i];

	return res;
}

// Multiply number by vector
template<class T>
inline NRvector<T> operator*(const T &a, const NRvector<T> &vec)
{
	return vec * a;
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != NULL) delete[] (v);
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_

template <class T>
class NRmatrix {
private:
	int nn;
	int mm;
	T **v;
public:
	NRmatrix();
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(int n, int m, const T &a);	//Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	NRmatrix & operator=(const NRmatrix &rhs);	//assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	~NRmatrix();
};

template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j,nel;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != NULL) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = nn>0 ? new T*[nn] : NULL;
			nel = mm*nn;
			if (v) v[0] = nel>0 ? new T[nel] : NULL;
			for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
		}
		for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* NRmatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline const T* NRmatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i];
}

template <class T>
inline int NRmatrix<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRmatrix<T>::ncols() const
{
	return mm;
}

template <class T>
void NRmatrix<T>::resize(int newn, int newm)
{
	int i,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
}

template <class T>
void NRmatrix<T>::assign(int newn, int newm, const T& a)
{
	int i,j,nel;
	if (newn != nn || newm != mm) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		v = nn>0 ? new T*[nn] : NULL;
		nel = mm*nn;
		if (v) v[0] = nel>0 ? new T[nel] : NULL;
		for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	}
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}


// basic type names (redefine if your bit lengths don't match)

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const NRvector<Int> vec_int_I;
typedef NRvector<Int> vec_int, vec_int_O, vec_int_IO;

typedef const NRvector<Uint> vec_uint_I;
typedef NRvector<Uint> vec_uint, vec_uint_O, vec_uint_IO;

typedef const NRvector<Llong> vec_llong_I;
typedef NRvector<Llong> vec_llong, vec_llong_O, vec_llong_IO;

typedef const NRvector<Ullong> vec_ullong_I;
typedef NRvector<Ullong> vec_ullong, vec_ullong_O, vec_ullong_IO;

typedef const NRvector<Char> vec_char_I;
typedef NRvector<Char> vec_char, vec_char_O, vec_char_IO;

typedef const NRvector<Char*> vec_charp_I;
typedef NRvector<Char*> vec_charp, vec_charp_O, vec_charp_IO;

typedef const NRvector<Uchar> vec_uchar_I;
typedef NRvector<Uchar> vec_uchar, vec_uchar_O, vec_uchar_IO;

typedef const NRvector<Doub> vec_doub_I;
typedef NRvector<Doub> vec_doub, vec_doub_O, vec_doub_IO;

typedef const NRvector<Doub*> vec_doubp_I;
typedef NRvector<Doub*> vec_doubp, vec_doubp_O, vec_doubp_IO;

typedef const NRvector<Complex> vec_complex_I;
typedef NRvector<Complex> vec_complex, vec_complex_O, vec_complex_IO;

typedef const NRvector<Bool> vec_bool_I;
typedef NRvector<Bool> vec_bool, vec_bool_O, vec_bool_IO;

// matrix types

typedef const NRmatrix<Int> MatInt_I;
typedef NRmatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const NRmatrix<Uint> MatUint_I;
typedef NRmatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const NRmatrix<Llong> MatLlong_I;
typedef NRmatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const NRmatrix<Ullong> MatUllong_I;
typedef NRmatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const NRmatrix<Char> MatChar_I;
typedef NRmatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const NRmatrix<Uchar> MatUchar_I;
typedef NRmatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const NRmatrix<Doub> MatDoub_I;
typedef NRmatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const NRmatrix<Bool> MatBool_I;
typedef NRmatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> Mat3DDoub_I;
typedef NRMat3d<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// Floating Point Exceptions for Microsoft compilers

#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions {
	turn_on_floating_exceptions() {
		int cw = _controlfp( 0, 0 );
		cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
		_controlfp( cw, MCW_EM );
	}
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */

}

#endif /* _NR3_H_ */

