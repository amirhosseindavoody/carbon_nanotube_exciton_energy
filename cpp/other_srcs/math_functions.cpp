#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <cstdlib> //defines EXIT_FAILURE and EXIT_SUCCESS

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "math_functions.h"

double math_functions::vector_norm(vector<double> &vec)
{
	double norm = 0;
	for(int i=0; i<vec.size(); i++)
	{
		norm += pow(vec[i], 2);
	}
	norm = sqrt(norm);
	return norm;
}


vector<double> math_functions::vector_sub(vector<double> &v1, vector<double> &v2)
{
	int len;
	if (v1.size() == v2.size())
	{
		len = v1.size();
	}
	else
	{
		cout << "unequal size of v1 and v2 in vector_sub" << endl;
		exit(EXIT_FAILURE);
	}
	vector<double> v_res(len);

	for (int i=0; i<len; i++)
	{
		v_res[i] = v1[i]-v2[i];
	}

	return v_res;
}

vector<double> math_functions::vector_add(vector<double> &v1, vector<double> &v2)
{
	int len;
	if (v1.size() == v2.size())
	{
		len = v1.size();
	}
	else
	{
		cout << "unequal size of v1 and v2 in vector_add" << endl;
		exit(EXIT_FAILURE);
	}
	vector<double> v_res(len);

	for (int i=0; i<len; i++)
	{
		v_res[i] = v1[i]+v2[i];
	}

	return v_res;
}

double math_functions::dot_product(vector<double> &v1, vector<double> &v2)
{
	int len;
	if (v1.size() == v2.size())
	{
		len = v1.size();
	}
	else
	{
		cout << "unequal size of v1 and v2 in dot_product" << endl;
		exit(EXIT_FAILURE);
	}
	
	double product = 0;

	for (int i=0; i<len; i++)
	{
		product += v1[i]*v2[i];
	}

	return product;
}

void math_functions::vector_clear(vector<int> &my_vector, int value)
{
	for (int i=0; i<my_vector.size(); i++)
	{
		my_vector[i] = value;
	}
}

void math_functions::vector_clear(vector<double> &my_vector, double value)
{
	for (int i=0; i<my_vector.size(); i++)
	{
		my_vector[i] = value;
	}
}

// Creates a vector with matlab style linspace numbering
vector<double> math_functions::linspace(double low, double high, int num)
{
	vector<double> my_vector(num);
	double step = (high - low) / static_cast<double>(num - 1);
	for (int i = 0; i < num; i++)
	{
		my_vector[i] = low;
		low += step;
	}
	return my_vector;
}


// Finds the index of the vector that has the number closest to but smaller than val.
// this function uses bisection method and relies on having an ascending or decending set of numbers in vec.
int math_functions::get_index(vector<double> &vec, double val)
{
	int ix_lower = 0;
	int ix_upper = vec.size()-1;
	int ix_mid;

	vector<double> tmp_vec(vec.size());
	for (int i=0; i<tmp_vec.size(); i++)
	{
		tmp_vec[i] = vec[i] - val;
	}

	if((tmp_vec[ix_lower]*tmp_vec[ix_upper]) >= 0)
	{
		if (tmp_vec[ix_lower] == 0) {return ix_lower;}
		else if(tmp_vec[ix_upper] == 0) {return ix_upper;}
		else
		{
			cout << "error in finding bisection root: value larger than bounds of the array" << endl;
			exit(EXIT_FAILURE);
		}
	}
	else if(tmp_vec[ix_lower] > 0)
	{
		while((ix_upper-ix_lower) > 1)
		{
			ix_mid = static_cast<int>(static_cast<double>(ix_upper + ix_lower) / 2.0);
			if(tmp_vec[ix_mid] > 0) {ix_lower = ix_mid;}
			else if(tmp_vec[ix_mid] < 0) {ix_upper = ix_mid;}
			else {return ix_mid;}
		}
		return ix_upper;
	}
	else
	{
		while((ix_upper-ix_lower) > 1)
		{
			ix_mid = static_cast<int>(static_cast<double>(ix_upper + ix_lower) / 2.0);
			if(tmp_vec[ix_mid] > 0) {ix_upper = ix_mid;}
			else if(tmp_vec[ix_mid] < 0) {ix_lower = ix_mid;}
			else {return ix_mid;}
		}
		return ix_lower;
	}
	
	cout << "error in finding bisection root!!!" << endl;
	exit(EXIT_FAILURE);
}