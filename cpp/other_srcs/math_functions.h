#ifndef math_functions_h
#define math_functions_h

#include <stdio.h>
#include <string>
#include <vector>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

using namespace std;

namespace math_functions
{
	double vector_norm(vector<double> &vec);
	vector<double> vector_sub(vector<double> &v1, vector<double> &v2);
	vector<double> vector_add(vector<double> &v1, vector<double> &v2);
	double dot_product(vector<double> &v1, vector<double> &v2);
	void vector_clear(vector<int> &my_vector, int value);
	void vector_clear(vector<double> &my_vector, double value);
	vector<double> linspace(double low, double high, int num);
	int get_index(vector<double> &vec, double val);
};

#endif //math_functions