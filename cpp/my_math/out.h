#ifndef _out_h_
#define _out_h_

#include "vector_matrix.h"

namespace nr
{
	// print complex matrix
	inline void print_mat(nr::mat_complex& X, std::string name="matrix_complex")
	{
		std::cout << name <<" = " << std::endl;
		for (int i=0; i<X.nrows(); i++)
		{
			for (int j=0; j<X.ncols(); j++)
			{
				std::cout << std::scientific;
				std::cout << std::showpos;
				std::cout << "\t" << X(i,j).real() << " " << X(i,j).imag() << " i";
			}
			std::cout << std::endl;
		}
	}

	// print double matrix
	inline void print_mat(nr::mat_doub& X, std::string name="matrix_double")
	{
		std::cout << name <<" = " << std::endl;
		for (int i=0; i<X.nrows(); i++)
		{
			for (int j=0; j<X.ncols(); j++)
			{
				std::cout << std::scientific;
				std::cout << std::showpos;
				std::cout << "\t" << X(i,j);
			}
			std::cout << std::endl;
		}
	}

	// print double vector
	inline void print_vec(nr::vec_doub& V, std::string name="vector_double")
	{
		std::cout << name <<" = " << std::endl;
		for (int i=0; i<V.size(); i++)
		{
			std::cout << std::scientific;
			std::cout << std::showpos;
			std::cout << "\t" << V(i) << std::endl;
		}
	}

}

#endif // _out_h_