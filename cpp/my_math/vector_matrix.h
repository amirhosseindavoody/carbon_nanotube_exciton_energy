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

// greatest common divisor
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


//*************************************************************************************************
//*************************************************************************************************
// Definition of vector class
//*************************************************************************************************
//*************************************************************************************************

#ifdef _USESTDVECTOR_
#define vector std::vector
#else

template<class T>
class matrix;

template <class T>
class vector {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	vector();
	explicit vector(int n);		// Zero-based array
	vector(int n, const T &a);	// Initialize to constant value
	vector(int n, const T *a);	// Initialize to array
	vector(const vector &rhs);	// Copy constructor
	vector & operator=(const vector &rhs);	//Copy assignment
	vector(vector &&rhs);	// Move constructor
	vector & operator=(vector &&rhs);	//Move assignment
	typedef T value_type; // make T available externally
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const; //i'th element
	inline T & operator()(const int i);	//i'th element
	inline const T & operator()(const int i) const; //i'th element
	inline int size() const; // gets size of the array
	inline const bool empty() const; // check if vector is empty
	void resize(int newn); // resize (contents not preserved)
	void assign(int newn, const T &a); // resize and assign a constant value
	vector operator+(const vector &other) const;	// addition operator
	vector operator-(const vector &other) const;	// subtract operator
	vector operator*(const T &a) const;	// multiply by number operator
	vector operator/(const T &a) const;	// divide by number operator
	T norm2() const; // norm of the vector
	template<class TT>
	friend vector<TT> operator*(const matrix<TT>& mat, const vector<TT>& vec); // matrix-vector multiplication
	template<class TT>
	friend TT dot_prod(const vector<TT> first, const vector<TT> second); // dot product of two vectors
	~vector();
};

// vector definitions

template <class T>
vector<T>::vector() : nn(0), v(NULL) {}

// Zero-based array
template <class T>
vector<T>::vector(int n) : nn(n), v(n>0 ? new T[n] : NULL) {}

// Initialize to constant value
template <class T>
vector<T>::vector(int n, const T& a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = a;
}

// Initialize to array
template <class T>
vector<T>::vector(int n, const T *a) : nn(n), v(n>0 ? new T[n] : NULL)
{
	for(int i=0; i<n; i++) v[i] = *a++;
}

// Copy constructor
template <class T>
vector<T>::vector(const vector<T> &rhs) : nn(rhs.nn), v(nn>0 ? new T[nn] : NULL)
{
	for(int i=0; i<nn; i++) v[i] = rhs[i];
}

// Move constructor
template <class T>
vector<T>::vector(vector<T> &&rhs) : nn(rhs.nn), v(nn>0 ? rhs.v : NULL)
{
	rhs.nn = 0;
	rhs.v = NULL;
}

// Copy assignment
template <class T>
vector<T> & vector<T>::operator=(const vector<T> &rhs)
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
vector<T> & vector<T>::operator=(vector<T> &&rhs)
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
inline T & vector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("vector subscript out of bounds");
}
#endif
	return v[i];
}

//i'th element
template <class T>
inline const T & vector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("vector subscript out of bounds");
}
#endif
	return v[i];
}

//i'th element
template <class T>
inline T & vector<T>::operator()(const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("vector subscript out of bounds");
}
#endif
	return v[i];
}

//i'th element
template <class T>
inline const T & vector<T>::operator()(const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("vector subscript out of bounds");
}
#endif
	return v[i];
}

// size of the vector
template <class T>
inline int vector<T>::size() const
{
	return nn;
}

// check if vector is empty
template <class T>
inline const bool vector<T>::empty() const
{
	return (v == NULL);
}

// resize the vector do not vector values
template <class T>
void vector<T>::resize(int newn)
{
	if (newn != nn) {
		if (v != NULL) delete[] (v);
		nn = newn;
		v = nn > 0 ? new T[nn] : NULL;
	}
}

// resize and assign a constant value
template <class T>
void vector<T>::assign(int newn, const T& a)
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
vector<T> vector<T>::operator+(const vector<T> &other) const
{
	if (v==NULL) throw("Trying to add empty vector!");
	if (other.v==NULL) throw("Trying to add empty vector!");
	if (nn != other.nn) throw("vector dimensions do not match for addition!");

	vector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i] + other[i];

	return res;
}

// subtract operator
template <class T>
vector<T> vector<T>::operator-(const vector<T> &other) const
{
	if (v==NULL) throw("error: empty vector!");
	if (other.v==NULL) throw("error: empty vector!");
	if (nn != other.nn) throw("vector dimensions mismatch!");

	vector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i] - other[i];

	return res;
}

// Multiply vector by number
template <class T>
vector<T> vector<T>::operator*(const T &a) const
{
	if (v==NULL) throw("Trying to multiply empty vector!");

	vector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = a * v[i];

	return res;
}

// Multiply number by vector
template<class T>
inline vector<T> operator*(const T &a, const vector<T> &vec)
{
	return vec * a;
}

// divide vector by number
template <class T>
vector<T> vector<T>::operator/(const T &a) const
{
	if (v==NULL) throw("error:empty vector!");

	vector<T> res(nn);
	for (int i = 0; i<nn; i++) res[i] = v[i]/a;

	return res;
}

// norm of the vector
template<class T>
T vector<T>::norm2() const
{
	if (v==NULL) throw("vector is empty!");

	T sum = 0;
	for (int i=0; i<nn; i++) sum += v[i]*v[i];

	return sqrt(sum);
}

// dot product of two vectors
template<class T>
inline T dot_prod(const vector<T> first, const vector<T> second)
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

// deconstructor
template <class T>
vector<T>::~vector()
{
	if (v != NULL) delete[] (v);
}

#endif //ifdef _USESTDVECTOR_

//*************************************************************************************************
//*************************************************************************************************
// definition of matrix class
//*************************************************************************************************
//*************************************************************************************************
template <class T>
class matrix {
private:
	int nn;
	int mm;
	T **v;
public:
	matrix(); // default constructor
	matrix(int n, int m);			// Zero-based array
	matrix(int n, int m, const T &a);	//Initialize to constant
	matrix(int n, int m, const T *a);	// Initialize to array
	matrix(const matrix &rhs);		// Copy constructor
	matrix(matrix &&rhs);		// move constructor
	matrix & operator=(const matrix &rhs);	//copy assignment
	matrix & operator=(matrix &&rhs);	//move assignment
	typedef T value_type; // make T available externally
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const; //subscripting: pointer to row i
	inline T& operator()(const int i, const int j);	//subscripting: reference to element at (i,j)
	inline const T& operator()(const int i, const int j) const; //subscripting: reference to element at (i,j)
	inline int dim1() const; // number of rows
	inline int dim2() const; // number of columns
	inline const bool empty() const; // check if matrix is empty
	void resize(int newn, int newm); // resize (contents not preserved)
	void assign(int newn, int newm, const T &a); // resize and assign a constant value
	matrix operator+(const matrix& other) const; // matrix-matrix addition
	matrix operator-(const matrix& other) const; // matrix-matrix addition
	matrix operator*(const T& num) const; // matrix-number multiplication
	matrix operator/(const T& num) const; // matrix-number multiplication
	matrix operator*(const matrix& other) const; // matrix-matrix multiplication
	template<class TT>
	friend vector<TT> operator*(const matrix<TT>& mat, const vector<TT>& vec); // matrix-vector multiplication
	~matrix();
};

// default constructor
template <class T>
matrix<T>::matrix() : nn(0), mm(0), v(NULL) {}

// Zero-based array
template <class T>
matrix<T>::matrix(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (int i=1;i<n;i++) v[i] = v[i-1] + m;
}

//Initialize to constant
template <class T>
matrix<T>::matrix(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

// Initialize to array
template <class T>
matrix<T>::matrix(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
	int i,j,nel=m*n;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< n; i++) v[i] = v[i-1] + m;
	for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

// Copy constructor
template <class T>
matrix<T>::matrix(const matrix<T> &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
	// std::cout << "matrix copy constructor!" << std::endl;
	int i,j,nel=mm*nn;
	if (v) v[0] = nel>0 ? new T[nel] : NULL;
	for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

// move constructor
template <class T>
matrix<T>::matrix(matrix<T> &&rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? rhs.v : NULL)
{
	// std::cout << "matrix move constructor!" << std::endl;
	rhs.nn = 0;
	rhs.mm = 0;
	rhs.v = NULL;
}

//copy assignment
template <class T>
matrix<T> & matrix<T>::operator=(const matrix<T> &rhs)
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
matrix<T> & matrix<T>::operator=(matrix<T> &&rhs)
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
inline T* matrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("matrix subscript out of bounds");
}
#endif
	return v[i];
}

//subscripting: pointer to row i
template <class T>
inline const T* matrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("matrix subscript out of bounds");
}
#endif
	return v[i];
}

//subscripting: reference to element at (i,j)
template <class T>
inline T& matrix<T>::operator()(const int i, const int j)
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("matrix subscript out of bounds");
}
if (j<0 || j>=mm) {
	throw("matrix subscript out of bounds");
}
#endif
	return v[i][j];
}

//subscripting: reference to element at (i,j)
template <class T>
inline const T& matrix<T>::operator()(const int i, const int j) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=nn) {
	throw("matrix subscript out of bounds");
}
if (j<0 || j>=mm) {
	throw("matrix subscript out of bounds");
}
#endif
	return v[i][j];
}

// number of rows
template <class T>
inline int matrix<T>::dim1() const
{
	return nn;
}

// number of columns
template <class T>
inline int matrix<T>::dim2() const
{
	return mm;
}

// check if matrix is empty
template <class T>
inline const bool matrix<T>::empty() const
{
	return (v == NULL);
} 

// resize (contents not preserved)
template <class T>
void matrix<T>::resize(int newn, int newm)
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
void matrix<T>::assign(int newn, int newm, const T& a)
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
matrix<T> matrix<T>::operator+(const matrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if ((nn != other.nn)||(mm!=other.mm)) throw("error: matrix dimensions mismatch!")
	matrix<T> res(nn,mm);

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
matrix<T> matrix<T>::operator-(const matrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if ((nn != other.nn)||(mm!=other.mm)) throw("error: matrix dimensions mismatch!")
	matrix<T> res(nn,mm);

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
matrix<T> matrix<T>::operator*(const T& num) const
{
	if (v == NULL) throw("error: empty matrix!");	
	matrix<T> res(nn,mm);

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
inline matrix<T> operator*(const T& num, const matrix<T>& mat)
{
	return mat*num;
}

// matrix-number division
template<class T>
matrix<T> matrix<T>::operator/(const T& num) const
{
	if (v == NULL) throw("error: empty matrix!");	
	matrix<T> res(nn,mm);

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
matrix<T> matrix<T>::operator*(const matrix<T>& other) const
{
	if (v == NULL) throw("error: empty matrix!");
	if (other.v == NULL) throw("error: empty matrix!");
	if (mm != other.nn) throw("multiplication error: size mismatch!");
	
	matrix<T> res(nn,other.mm);
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
matrix<T>::~matrix()
{
	if (v != NULL) {
		delete[] (v[0]);
		delete[] (v);
	}
}

// matrix-vector multiplication
template<class T>
vector<T> operator*(const matrix<T>& mat, const vector<T>& vec)
{
	#ifdef _CHECKEMPTY_
	if (vec.empty()) throw("error: empty vector!");
	if (mat.empty()) throw("error: empty matrix!");
	#endif // end _CHECKEMPTY_

	#ifdef _CHECKDIMENSIONS_
	if (vec.size() != mat.dim2()) throw("multiplication error: size mismatch!");
	#endif // end _CHECKDIMENSIONS_
	
	vector<T> res(mat.dim1());
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

//*************************************************************************************************
//*************************************************************************************************
// definition of 3D matrix class
//*************************************************************************************************
//*************************************************************************************************
template <class T>
class matrix3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	matrix3d(); // constructor
	matrix3d(const int n, const int m, const int k); // zero-based array
	matrix3d(const int n, const int m, const int k, const T &a ); // initiallize to constant value
	matrix3d(const matrix3d &rhs); // copy constructor
	matrix3d(matrix3d &&rhs); // move constructor
	inline matrix3d & operator=(const matrix3d &rhs);	//copy assignment
	inline matrix3d & operator=(matrix3d &&rhs);	//move assignment
	inline T** operator()(const int i);	//subscripting: pointer to row i
	inline T& operator()(const int i, const int j, const int k);	//subscripting: reference to element at (i,j,k)
	inline const T& operator()(const int i, const int j, const int k) const; //subscripting: reference to element at (i,j,k)
	inline const T* const * operator()(const int i) const; //subscripting: pointer to row i
	inline int dim1() const; //size of dimension 1
	inline int dim2() const; //size of dimension 2
	inline int dim3() const; //size of dimension 3
	void resize(const int newn, const int newm, const int newk); // resize (contents not preserved)
	void assign(const int newn, const int newm, const int newk, const T &a); // resize and assign a constant value
	~matrix3d(); // deconstructor
};

// constructor
template <class T>
matrix3d<T>::matrix3d(): nn(0), mm(0), kk(0), v(NULL) {}

// zero-based array
template <class T>
matrix3d<T>::matrix3d(const int n, const int m, const int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(int j=1; j<m; j++)
	{
		v[0][j] = v[0][j-1] + k;
	}
	for(int i=1; i<n; i++)
	{
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(int j=1; j<m; j++)
		{
			v[i][j] = v[i][j-1] + k;
		}
	}
}

// initiallize to constant value
template <class T>
matrix3d<T>::matrix3d(const int n, const int m, const int k, const T &a) : nn(n), mm(m), kk(k), v(new T**[n])
{
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(int j=1; j<m; j++)
	{
		v[0][j] = v[0][j-1] + k;
	}
	for(int i=1; i<n; i++)
	{
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(int j=1; j<m; j++)
		{
			v[i][j] = v[i][j-1] + k;
		}
	}
	for (int i=0; i<nn; i++)
	{
		for (int j=0; j<mm; j++)
		{
			for (int l=0; l<kk; l++)
			{
				v[i][j][l] = a;
			}
		}
	}
}

// copy constructor
template<class T>
matrix3d<T>::matrix3d(const matrix3d<T> &rhs): nn(rhs.nn), mm(rhs.mm), kk(rhs.kk), v(nn>0 ? new T**[nn] : NULL)
{
	if (v)
	{
		v[0] = nn*mm>0 ? new T*[nn*mm] : NULL;
	}
	if (v[0])
	{
		v[0][0] = nn*mm*kk>0 ? new T[nn*mm*kk] : NULL;
	}
	for(int j=1; j<mm; j++)
	{
		v[0][j] = v[0][j-1] + kk;
	}
	for (int i=1; i< nn; i++)
	{
		v[i] = v[i-1] + mm;
		v[i][0] = v[i-1][0] + mm*kk;
		for (int j=1; j<mm; j++)
		{
			v[i][j] = v[i][j-1]+kk;
		}
	}
	for (int i=0; i<nn; i++)
	{
		for (int j=0; j<mm; j++)
		{
			for (int k=0; k<kk; k++)
			{
				v[i][j][k] = rhs(i,j,k);
			}
		}
	}
}

// move constructor
template<class T>
matrix3d<T>::matrix3d(matrix3d<T> &&rhs): nn(rhs.nn), mm(rhs.mm), kk(rhs.kk), v(nn>0 ? rhs.v : NULL)
{
	rhs.nn = 0;
	rhs.mm = 0;
	rhs.kk = 0;
	rhs.v = NULL;
}

// copy assignment
template <class T>
inline matrix3d<T> & matrix3d<T>::operator=(const matrix3d<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.dim1() || mm != rhs.dim2() || kk!=rhs.dim3())
		{
			if (v != NULL)
			{
				delete[] (v[0][0]);
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.dim1();
			mm=rhs.dim2();
			kk=rhs.dim3();
			v = nn>0 ? new T**[nn] : NULL;
			if (v)
			{
				v[0] = nn*mm>0 ? new T*[nn*mm] : NULL;
			}
			if (v[0])
			{
				v[0][0] = nn*mm*kk>0 ? new T[nn*mm*kk] : NULL;
			}
			for(int j=1; j<mm; j++)
			{
				v[0][j] = v[0][j-1] + kk;
			}
			for (int i=1; i< nn; i++)
			{
				v[i] = v[i-1] + mm;
				v[i][0] = v[i-1][0] + mm*kk;
				for (int j=1; j<mm; j++)
				{
					v[i][j] = v[i][j-1]+kk;
				}
			}
		}

		for (int i=0; i<nn; i++)
		{
			for (int j=0; j<mm; j++)
			{
				for (int k=0; k<kk; k++)
				{
					v[i][j][k] = rhs(i,j,k);
				}
			}
		}
	}
	return *this;
}

// move assignment
template <class T>
inline matrix3d<T> & matrix3d<T>::operator=(matrix3d<T> &&rhs)
{

	if (this != &rhs)
	{
		if (v != NULL)
		{
			delete[] (v[0][0]);
			delete[] (v[0]);
			delete[] (v);
		}
		nn=rhs.dim1();
		mm=rhs.dim2();
		kk=rhs.dim3();
		if (rhs.v)
		{
			v = rhs.v;
		}
		else
		{
			v = NULL;
		}

		rhs.nn = 0;
		rhs.mm = 0;
		rhs.kk = 0;
		rhs.v = NULL;
	}

	return *this;
}

//subscripting: pointer to row i
template <class T>
inline T** matrix3d<T>::operator()(const int i) 
{
	return v[i];
}

//subscripting: pointer to row i
template <class T>
inline const T* const * matrix3d<T>::operator()(const int i) const
{
	return v[i];
}

//subscripting: reference to element at (i,j,k)
template <class T>
inline T& matrix3d<T>::operator()(const int i, const int j, const int k)
{
	#ifdef _CHECKBOUNDS_
	if (i<0 || i>=dim1())
	{
		throw("matrix subscript out of bounds");
	}
	if (j<0 || j>=dim2())
	{
		throw("matrix subscript out of bounds");
	}
	if (k<0 || k>=dim3())
	{
		throw("matrix subscript out of bounds");
	}
	#endif
	return v[i][j][k];
}

//subscripting: reference to element at (i,j,k)
template <class T>
inline const T& matrix3d<T>::operator()(const int i, const int j, const int k) const
{
	#ifdef _CHECKBOUNDS_
	if (i<0 || i>=dim1())
	{
		throw("matrix subscript out of bounds");
	}
	if (j<0 || j>=dim2())
	{
		throw("matrix subscript out of bounds");
	}
	if (k<0 || k>=dim3())
	{
		throw("matrix subscript out of bounds");
	}
	#endif
	return v[i][j][k];
}

//size of dimension 1
template <class T>
inline int matrix3d<T>::dim1() const
{
	return nn;
}

//size of dimension 2
template <class T>
inline int matrix3d<T>::dim2() const
{
	return mm;
}

//size of dimension 3
template <class T>
inline int matrix3d<T>::dim3() const
{
	return kk;
}

// resize (contents not preserved)
template <class T>
void matrix3d<T>::resize(const int newn, const int newm, const int newk)
{
	if (newn != nn || newm != mm || newk != kk)
	{
		if (v != NULL) {
			delete[] (v[0][0]);
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		kk = newk;
		v = nn>0 ? new T**[nn] : NULL;
		if (v)
		{
			v[0] = nn*mm>0 ? new T*[nn*mm] : NULL;
		}
		if (v[0])
		{
			v[0][0] = nn*mm*kk>0 ? new T[nn*mm*kk] : NULL;
		}
		for(int j=1; j<mm; j++)
		{
			v[0][j] = v[0][j-1] + kk;
		}
		for (int i=1; i< nn; i++)
		{
			v[i] = v[i-1] + mm;
			v[i][0] = v[i-1][0] + mm*kk;
			for (int j=1; j<mm; j++)
			{
				v[i][j] = v[i][j-1]+kk;
			}
		}
	}
}

// resize and assign a constant value
template <class T>
void matrix3d<T>::assign(const int newn, const int newm, const int newk, const T& a)
{
	if (newn != nn || newm != mm || newk != kk)
	{

		if (v != NULL) {
			delete[] (v[0][0]);
			delete[] (v[0]);
			delete[] (v);
		}
		nn = newn;
		mm = newm;
		kk = newk;
		v = nn>0 ? new T**[nn] : NULL;
		if (v)
		{
			v[0] = nn*mm>0 ? new T*[nn*mm] : NULL;
		}
		if (v[0])
		{
			v[0][0] = nn*mm*kk>0 ? new T[nn*mm*kk] : NULL;
		}
		for(int j=1; j<mm; j++)
		{
			v[0][j] = v[0][j-1] + kk;
		}
		for (int i=1; i< nn; i++)
		{
			v[i] = v[i-1] + mm;
			v[i][0] = v[i-1][0] + mm*kk;
			for (int j=1; j<mm; j++)
			{
				v[i][j] = v[i][j-1]+kk;
			}
		}
	}

	for (int i=0; i< nn; i++)
	{
		for (int j=0; j<mm; j++)
		{
			for (int k=0; k<kk; k++)
			{
				v[i][j][k] = a;
			}
		}
	}
}

// deconstrocture
template <class T>
matrix3d<T>::~matrix3d()
{
	if (v != NULL) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

//*************************************************************************************************
//*************************************************************************************************
// vector types
//*************************************************************************************************
//*************************************************************************************************

typedef const vector<Int> vec_int_I;
typedef vector<Int> vec_int, vec_int_O, vec_int_IO;

typedef const vector<Uint> vec_uint_I;
typedef vector<Uint> vec_uint, vec_uint_O, vec_uint_IO;

typedef const vector<Llong> vec_llong_I;
typedef vector<Llong> vec_llong, vec_llong_O, vec_llong_IO;

typedef const vector<Ullong> vec_ullong_I;
typedef vector<Ullong> vec_ullong, vec_ullong_O, vec_ullong_IO;

typedef const vector<Char> vec_char_I;
typedef vector<Char> vec_char, vec_char_O, vec_char_IO;

typedef const vector<Char*> vec_charp_I;
typedef vector<Char*> vec_charp, vec_charp_O, vec_charp_IO;

typedef const vector<Uchar> vec_uchar_I;
typedef vector<Uchar> vec_uchar, vec_uchar_O, vec_uchar_IO;

typedef const vector<Doub> vec_doub_I;
typedef vector<Doub> vec_doub, vec_doub_O, vec_doub_IO;

typedef const vector<Doub*> vec_doubp_I;
typedef vector<Doub*> vec_doubp, vec_doubp_O, vec_doubp_IO;

typedef const vector<cmplx> vec_complex_I;
typedef vector<cmplx> vec_complex, vec_complex_O, vec_complex_IO;

typedef const vector<Bool> vec_bool_I;
typedef vector<Bool> vec_bool, vec_bool_O, vec_bool_IO;

//*************************************************************************************************
//*************************************************************************************************
// matrix types
//*************************************************************************************************
//*************************************************************************************************

typedef const matrix<Int> mat_int_I;
typedef matrix<Int> mat_int, mat_int_O, mat_int_IO;

typedef const matrix<Uint> mat_uint_I;
typedef matrix<Uint> mat_uint, mat_uint_O, mat_uint_IO;

typedef const matrix<Llong> mat_llong_I;
typedef matrix<Llong> mat_llong, mat_llong_O, mat_llong_IO;

typedef const matrix<Ullong> mat_ullong_I;
typedef matrix<Ullong> mat_ullong, mat_ullong_O, mat_ullong_IO;

typedef const matrix<Char> mat_char_I;
typedef matrix<Char> mat_char, mat_char_O, mat_char_IO;

typedef const matrix<Uchar> mat_uchar_I;
typedef matrix<Uchar> mat_uchar, mat_uchar_O, mat_uchar_IO;

typedef const matrix<Doub> mat_doub_I;
typedef matrix<Doub> mat_doub, mat_doub_O, mat_doub_IO;

typedef const matrix<cmplx> mat_complex_I;
typedef matrix<cmplx> mat_complex, mat_complex_O, mat_complex_IO;

typedef const matrix<Bool> mat_bool_I;
typedef matrix<Bool> mat_bool, mat_bool_O, mat_bool_IO;

//*************************************************************************************************
//*************************************************************************************************
// 3D matrix types
//*************************************************************************************************
//*************************************************************************************************

typedef const matrix3d<Doub> mat3d_doub_I;
typedef matrix3d<Doub> mat_3d_doub, mat3d_doub_O, mat3d_doub_IO;

typedef const matrix3d<cmplx> mat3d_complex_I;
typedef matrix3d<cmplx> mat3d_complex, mat3d_complex_O, mat3d_complex_IO;

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

