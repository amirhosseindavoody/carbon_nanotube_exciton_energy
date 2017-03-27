#ifndef d_vector_h
#define d_vector_h

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <valarray>

#include <Eigen/Dense>

namespace my_math
{
	class c_matrix;
	class c_vector;
	class d_matrix;

	// complex vector class
	class d_vector
	{
	private:
		// define element type
		typedef double element_t;

		int row;
		const int column=1;
		Eigen::Matrix<element_t, Eigen::Dynamic, 1> mat;

	public:
		// constructor
		d_vector(int n)
		{
			row = n;
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
		element_t& operator()(int n)
		{
			assert(n >= 0 && n < row);
			return mat(n,0);
		}
		
		// get elements by the index
		const element_t& operator()(int n) const
		{
			assert(n >= 0 && n < row);
			return mat(n,0);
		}
	};

};

#endif //d_vector