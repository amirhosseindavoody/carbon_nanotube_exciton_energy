#include <stdio.h>
#include <iostream>
#include <string>
#include <cstdlib> //defines EXIT_FAILURE and EXIT_SUCCESS

#include "write_log.h"
#include "math_functions.h"
#include "tableElem.h"

using namespace math_functions;

// Creates table element object
tableElem::tableElem(double distance, double angle, double forster_const, int tube_idx, int seg_idx)
{
	r = distance;
	theta = angle;
	gamma = forster_const;
	tubeidx = tube_idx;
	segidx = seg_idx;

	set_rate();
}

void tableElem::set_rate()
{
	rate = abs(gamma*cos(theta) / pow(r, 6));
}


// Gets the total transition rate based on gamma, r, and theta
double tableElem::getRate()
{
	return rate;
}


// Gets r value
double tableElem::getr()
{
	return r;
}


// Gets theta value
double tableElem::getTheta()
{
	return theta;
}

// Gets tube number
int tableElem::getTubeidx()
{
	return tubeidx;
}


// Gets the segment number
int tableElem::getSegidx()
{
	return segidx;
}


