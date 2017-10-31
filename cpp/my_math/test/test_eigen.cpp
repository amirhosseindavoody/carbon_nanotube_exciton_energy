#include <iostream>
#include <math.h>
#include <complex>

#include "nr.h"

int main(int argc, char *argv[])
{

	std::cout << std::endl << "testing eigen value solver for complex hermitian matrix:" << std::endl;

	int n = 2;
	nr::mat_complex X(n,n,0.0);
	nr::mat_complex eigvec(n,n,0.0);
	nr::vec_doub eigval(n,0.0);

	X(0,0) = std::complex<double>(+1.0,+0.0);
	X(0,1) = std::complex<double>(+0.0,+2.0);
	X(1,0) = std::complex<double>(+0.0,-2.0);
	X(1,1) = std::complex<double>(-3.0,+0.0);
	nr::print(X, "X");

	nr::eig_sym(eigval, eigvec, X);
	
	// change phase of eigen vectors for better appearance
	for (int i=0; i<eigvec.nrows(); i++)
	{
		for (int j=0; j<eigvec.ncols(); j++)
		{
			eigvec(i,j) = eigvec(i,j)/eigvec(eigvec.nrows()-1,j);
		}
	}

	// check eigen vector normalization and orthogonality
	for (int i=0; i<eigvec.ncols(); i++)
	{
		std::complex<double> sum(0.0,0.0);
		for (int j=0; j<eigvec.nrows(); j++)
		{
			sum += std::conj(eigvec(j,i)) * eigvec(j,0);
		}
		std::cout << i << " , sum = " << std::abs(sum) << std::endl;
	}

	nr::print(eigval, "eigval");
	nr::print(eigvec, "eigvec");

	nr::mat_complex lambda(n,n,0.0);
	for (int i=0; i<lambda.nrows(); i++)
	{
		lambda(i,i) = eigval(i);
	}
	nr::print(lambda, "lambda");

	nr::mat_complex test(X.nrows(),X.ncols());
	test = (X*eigvec-eigvec*lambda);
	nr::print(test, "test");

	std::cout << "if the matrix test is approximately zero then the eigen value solution was successful!" << std::endl;
	
	return 0;
}