#ifndef my_math_h
#define my_math_h

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <valarray>

#include <Eigen/Dense>

#include "c_vector.h"
#include "c_matrix.h"
#include "d_matrix.h"
#include "d_vector.h"

#include "multiplication.h"

namespace my_math
{
	int gcd(int a, int b);
	
	template<class myType>
	myType my_norm(std::valarray<myType> a)
	{
		return sqrt((a*a).sum());
	}

};

#endif //my_math