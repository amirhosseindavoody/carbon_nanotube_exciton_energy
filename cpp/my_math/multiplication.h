#ifndef multiplication_h
#define multiplication_h

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <valarray>

#include <Eigen/Dense>

#include "my_math.h"

#include "c_vector.h"
#include "c_matrix.h"
#include "d_matrix.h"
#include "d_vector.h"

namespace my_math
{

	// matrix * vector multiplication
	c_vector& matmul(const c_matrix& A, const c_vector& x)
	{
		typedef std::complex<double> element_t;

		assert(A.size2() == x.size1());
		c_vector result(A.size1());
		result.mat = A.mat * x.mat;
		return result;
	};

};

#endif //multiplication_h