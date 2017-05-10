#ifndef _vector_matrix_h_
#define _vector_matrix_h_

// check bounds when accessing elements in vectors and matrices
#define _CHECKBOUNDS_ 1	
// check if vectors and matrices are empty when operating on these objects
#define _CHECKEMPTY_ 1
// check if dimensions matches in algebraic operations on matrices and vectors
#define _CHECKDIMENSIONS_ 1 
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

// // greatest common divisor
inline int gcd ( const int a, const int b ) 
{ 
  int c;
  int aa = a;
  int bb = b;
  while ( aa != 0 ) 
  { 
     c = aa; 
     aa = bb%aa; 
     bb = c; 
  } 
  return bb; 
}


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

template<class T>
class NRmatrix;

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
	inline const T & operator[](const int i) const; //i'th element
	inline T & operator()(const int i);	//i'th element
	inline const T & operator()(const int i) const; //i'th element
	inline int size() const; // gets size of the array
	inline const bool empty() const; // check if vector is empty
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	NRvector operator+(const NRvector &other) const;	// addition operator
	NRvector operator-(const NRvector &other) const;	// subtract operator
	NRvector operator*(const T &a) const;	// multiply by number operator
	NRvector operator/(const T &a) const;	// divide by number operator
	T norm2() const; // norm of the vector
	template<class TT>
	friend NRvector<TT> operator*(const NRmatrix<TT>& mat, const NRvector<TT>& vec); // matrix-vector multiplication
	template<class TT>
	friend TT dot_prod(const NRvector<TT> first, const NRvector<TT> second); // dot product of two vectors
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
{
	// std::cout << "vector move assignment!" << std::endl;
	if (this != &rhs)
	{
		if (v != NULL) delete [] (v);
		nn=rhs.nn;
		v= nn>0 ? rhs.v : NULL;
		rhs.nn = 0;
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

//i'th element
template <class T>
inline T & NRvector<T>::operator()(const int i)	//subscripting
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
inline const T & NRvector<T>::operator()(const int i) const	//subscripting
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

// check if vector is empty
template <class T>
inline const bool NRvector<T>::empty() const
{
	return (v == NULL);
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

// addition operator
template <class T>
NRvector<T> NRvector<T>::operator+(const NRvector<T> &other) const
{
	if (v==NULL) throw("Trying to add empty vector!");
	if (other.v==NULL) throw("Trying to add empty vector!");
	if (nn != other.nn) throw("NRvector dimensions do not match for addition!");

	NRvector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i] + other[i];

	return res;
}

// subtract operator
template <class T>
NRvector<T> NRvector<T>::operator-(const NRvector<T> &other) const
{
	if (v==NULL) throw("error: empty vector!");
	if (other.v==NULL) throw("error: empty vector!");
	if (nn != other.nn) throw("vector dimensions mismatch!");

	NRvector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i] - other[i];

	return res;
}

// Multiply vector by number
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

// divide vector by number
template <class T>
NRvector<T> NRvector<T>::operator/(const T &a) const
{
	if (v==NULL) throw("error:empty vector!");

	NRvector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i]/a;

	return res;
}

// norm of the vector
template<class T>
T NRvector<T>::norm2() const
{
	if (v==NULL) throw("vector is empty!");

	T sum = 0;
	for (int i=0; i<nn; i++) sum += v[i]*v[i];

	return sqrt(sum);
}

// dot product of two vectors
template<class T>
inline T dot_prod(const NRvector<T> first, const NRvector<T> second)
{
	if (first.v == NULL) throw("error: vector is empty!");
	if (second.v == NULL) throw("error: vector is empty!");
	if (first.nn != second.nn) throw("error: vector dimensions mismatch!")

	T sum = static_cast<T>(0);
	for (int i=0; i<first.nn; i++)
	{
		sum += first.v[i]*second.v[i];
	}
	return sum;
}

template <class T>
NRvector<T>::~NRvector()
{
	if (v != NULL) delete[] (v);
}

// end of NRvector definitions

#endif //ifdef _USESTDVECTOR_


// definition of matrix class
template <class T>
class NRmatrix {
private:
	int nn;
	int mm;
	T **v;
public:
	NRmatrix(); // default constructor
	NRmatrix(int n, int m);			// Zero-based array
	NRmatrix(int n, int m, const T &a);	//Initialize to constant
	NRmatrix(int n, int m, const T *a);	// Initialize to array
	NRmatrix(const NRmatrix &rhs);		// Copy constructor
	NRmatrix(NRmatrix &&rhs);		// move constructor
	NRmatrix & operator=(const NRmatrix &rhs);	//copy assignment
	NRmatrix & operator=(NRmatrix &&rhs);	//move assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const; //subscripting: pointer to row i
	inline T& operator()(const int i, const int j);	//subscripting: reference to element at (i,j)
	inline const T& operator()(const int i, const int j) const; //subscripting: reference to element at (i,j)
	inline int nrows() const; // number of rows
	inline int ncols() const; // number of columns
	inline const bool empty() const; // check if matrix is empty
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	NRmatrix operator+(const NRmatrix& other) const; // matrix-matrix addition
	NRmatrix operator-(const NRmatrix& other) const; // matrix-matrix addition
	NRmatrix operator*(const T& num) const; // matrix-number multiplication
	NRmatrix operator/(const T& num) const; // matrix-number multiplication
	NRmatrix operator*(const NRmatrix& other) const; // matrix-matrix multiplication
	template<class TT>
	friend NRvector<TT> operator*(const NRmatrix<TT>& mat, const NRvector<TT>& vec); // matrix-vector multiplication
	~NRmatrix();
};

// default constructor
template <class T>
NRmatrix<T>::NRmatrix() : nn(0), mm(0), v(NULL) {}

// Zero-based array
template <class T>
NRmatrix<T>::NRmatrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (int i=1;i<n;i++) v[i] = v[i-1] + m;
}

//Initialize to constant
template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

// Initialize to array
template <class T>
NRmatrix<T>::NRmatrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

// Copy constructor
template <class T>
NRmatrix<T>::NRmatrix(const NRmatrix<T> &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	// std::cout << "matrix copy constructor!" << std::endl;
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

// move constructor
template <class T>
NRmatrix<T>::NRmatrix(NRmatrix<T> &&rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? rhs.v : NULL)
{
	// std::cout << "matrix move constructor!" << std::endl;
	rhs.nn = 0;
	rhs.mm = 0;
	rhs.v = NULL;
}

//copy assignment
template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(const NRmatrix<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	// std::cout << "matrix copy assignment!" << std::endl;
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

//move assignment
template <class T>
NRmatrix<T> & NRmatrix<T>::operator=(NRmatrix<T> &&rhs)
{
	// std::cout << "matrix move assignment!" << std::endl;
	if (this != &rhs) {
		if (v != NULL) {
			delete[] (v[0]);
			delete[] (v);
		}
		nn=rhs.nn;
		mm=rhs.mm;
		if (rhs.v) v = rhs.v;
		else v = NULL;

		rhs.nn = 0;
		rhs.mm = 0;
		rhs.v = NULL;
	}
	return *this;
}

//subscripting: pointer to row i
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

//subscripting: pointer to row i
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

//subscripting: reference to element at (i,j)
template <class T>
inline T& NRmatrix<T>::operator()(const int i, const int j)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
if (j<0 || j>=mm) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i][j];
}

//subscripting: reference to element at (i,j)
template <class T>
inline const T& NRmatrix<T>::operator()(const int i, const int j) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("NRmatrix subscript out of bounds");
}
if (j<0 || j>=mm) {
	throw("NRmatrix subscript out of bounds");
}
#endif
	return v[i][j];
}

// number of rows
template <class T>
inline int NRmatrix<T>::nrows() const
{
	return nn;
}

// number of columns
template <class T>
inline int NRmatrix<T>::ncols() const
{
	return mm;
}

// check if matrix is empty
template <class T>
inline const bool NRmatrix<T>::empty() const
{
	return (v == NULL);
} 

// resize (contents not preserved)
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

// resize and assign a constant value
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

// matrix-matrix addition
template<class T>
NRmatrix<T> NRmatrix<T>::operator+(const NRmatrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if ((nn != other.nn)||(mm!=other.mm)) throw("error: matrix dimensions mismatch!")
	NRmatrix<T> res(nn,mm);

	for (int i=0; i<res.nn; i++)
	{
		for (int j=0; j<res.mm; j++)
		{
			res.v[i][j] = v[i][j]+other.v[i][j];
		}	
	}
	return res;
}

// matrix-matrix subtraction
template<class T>
NRmatrix<T> NRmatrix<T>::operator-(const NRmatrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if ((nn != other.nn)||(mm!=other.mm)) throw("error: matrix dimensions mismatch!")
	NRmatrix<T> res(nn,mm);

	for (int i=0; i<res.nn; i++)
	{
		for (int j=0; j<res.mm; j++)
		{
			res.v[i][j] = v[i][j]-other.v[i][j];
		}	
	}
	return res;
}

// matrix-number multiplication
template<class T>
NRmatrix<T> NRmatrix<T>::operator*(const T& num) const
{
	if (v == NULL) throw("error: empty matrix!");	
	NRmatrix<T> res(nn,mm);

	for (int i=0; i<res.nn; i++)
	{
		for (int j=0; j<res.mm; j++)
		{
			res.v[i][j] = num*v[i][j];
		}	
	}
	return res;
}

template<class T>
inline NRmatrix<T> operator*(const T& num, const NRmatrix<T>& mat)
{
	return mat*num;
}

// matrix-number division
template<class T>
NRmatrix<T> NRmatrix<T>::operator/(const T& num) const
{
	if (v == NULL) throw("error: empty matrix!");	
	NRmatrix<T> res(nn,mm);

	for (int i=0; i<res.nn; i++)
	{
		for (int j=0; j<res.mm; j++)
		{
			res.v[i][j] = v[i][j]/num;
		}	
	}
	return res;
}

// matrix-matrix multiplication
template<class T>
NRmatrix<T> NRmatrix<T>::operator*(const NRmatrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if (mm != other.nn) throw("multiplication error: size mismatch!");
	
	NRmatrix<T> res(nn,other.mm);
	T sum = static_cast<T>(0);

	for (int i=0; i<res.nn; i++)
	{
		for (int j=0; j<res.mm; j++)
		{
			sum = static_cast<T>(0);
			for (int k=0; k<mm; k++)
			{
				sum += v[i][k]*other.v[k][j];
			}
			res.v[i][j] = sum;
		}	
	}
	return res;
}

// deconstructor
template <class T>
NRmatrix<T>::~NRmatrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

// matrix-vector multiplication
template<class T>
NRvector<T> operator*(const NRmatrix<T>& mat, const NRvector<T>& vec)
{
	#ifdef _CHECKEMPTY_
	if (vec.empty()) throw("error: empty vector!");
	if (mat.empty()) throw("error: empty matrix!");
	#endif // end _CHECKEMPTY_

	#ifdef _CHECKDIMENSIONS_
	if (vec.size() != mat.ncols()) throw("multiplication error: size mismatch!");
	#endif // end _CHECKDIMENSIONS_
	
	NRvector<T> res(mat.nrows());
	T sum = static_cast<T>(0);

	for (int i=0; i<res.size(); i++)
	{
		sum = static_cast<T>(0);
		for (int j=0; j<vec.size(); j++)
		{
			sum += mat(i,j)*vec(j);
		}
		res(i) = sum;
	}
	return res;
}

// definition of 3D matrix class
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

typedef const NRmatrix<Int> mat_int_I;
typedef NRmatrix<Int> mat_int, mat_int_O, mat_int_IO;

typedef const NRmatrix<Uint> mat_uint_I;
typedef NRmatrix<Uint> mat_uint, mat_uint_O, mat_uint_IO;

typedef const NRmatrix<Llong> mat_llong_I;
typedef NRmatrix<Llong> mat_llong, mat_llong_O, mat_llong_IO;

typedef const NRmatrix<Ullong> mat_ullong_I;
typedef NRmatrix<Ullong> mat_ullong, mat_ullong_O, mat_ullong_IO;

typedef const NRmatrix<Char> mat_char_I;
typedef NRmatrix<Char> mat_char, mat_char_O, mat_char_IO;

typedef const NRmatrix<Uchar> mat_uchar_I;
typedef NRmatrix<Uchar> mat_uchar, mat_uchar_O, mat_uchar_IO;

typedef const NRmatrix<Doub> mat_doub_I;
typedef NRmatrix<Doub> mat_doub, mat_doub_O, mat_doub_IO;

typedef const NRmatrix<Complex> mat_complex_I;
typedef NRmatrix<Complex> mat_complex, mat_complex_O, mat_complex_IO;

typedef const NRmatrix<Bool> mat_bool_I;
typedef NRmatrix<Bool> mat_bool, mat_bool_O, mat_bool_IO;

// 3D matrix types

typedef const NRMat3d<Doub> mat_3d_doub_I;
typedef NRMat3d<Doub> mat_3d_doub, mat_3d_doub_O, mat_3d_doub_IO;

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

#endif // _vector_matrix_h_

