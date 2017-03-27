#ifndef d_matrix_h
#define d_matrix_h

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <valarray>

#include <Eigen/Dense>

namespace my_math
{
	class c_matrix;
	class c_vector;
	class d_vector;

// complex matrix class
	class d_matrix
	{
	private:
		// define element type
		typedef double element_t;

		int row, column;
		Eigen::Matrix<element_t, Eigen::Dynamic, Eigen::Dynamic> mat;

	public:
		// constructor
		d_matrix(int n, int m)
		{
			row = n;
			column = m;
			mat.resize(row,column);
			zero();
		}

		// zeros all the elements
		void zero()
		{
		for (int i=0; i<row; i++)
			for (int j=0; j<column; j++)
				mat(i,j) = 0.e0;
		}

		// gets the number of rows
		int size1() const
		{
			return row;
		}

		// gets the number of columns
		int size2() const
		{
			return column;
		}

		// get elements by the index
		element_t& operator()(int n, int m)
		{
			assert(n >= 0 && n < row);
			assert(m >= 0 && m < column);
			return mat(n,m);
		}

		// get elements by the index
		const element_t& operator()(int n, int m) const
		{
			assert(n >= 0 && n < row);
			assert(m >= 0 && m < column);
			return mat(n,m);
		}	
	};

};

#endif //d_matrix