#include <stdio.h>
#include <iostream>

#include "segment.h"
#include "math_functions.h"

using namespace math_functions;

// Determines if the segment has an exciton of the same type as the passed exciton
segment::segment(int number, vector<double> first_point, vector<double> second_point)
{
	segment_number = number;
	
	point1 = first_point;
	point2 = second_point;
	
	point_m = vector<double>(3);
	for (int i=0; i<3; i++)
	{
		point_m[i] = (point1[i]+point2[i])/2.0;
	}

	tbl = vector<tableElem>(0);
	rateVec = vector<double>(0);
}

bool segment::hasExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		exit(EXIT_FAILURE);
	}

	//if energy is 1 and that slot is occupied return true
	if (currEn == 1 && !(ex1 == nullptr))
	{
		return true;
	}
	//if energy is 2 and that slot is occupied return true
	if (currEn == 2 && !(ex2 == nullptr))
	{
		return true;
	}
	return false; //slot is empty
}

// Sets the exciton in the correct slot to the exciton passed to function.
bool segment::setExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	if (!this->hasExciton(e))
	{
		if (currEn == 1)
		{
			ex1 = e;
			return true;
		}
		if (currEn == 2)
		{
			ex2 = e;
			return true;
		}
	}
	return false;
}

// Removes the exciton of the correct type
bool segment::removeExciton(shared_ptr<exciton> e)
{
	int currEn = e->getEnergy();
	//parameter check
	if (currEn != 1 && currEn != 2)
	{
		cout << "Exciton passed to \"hasExciton\" has not been initialized correctly.\n";
		system("pause");
		exit(EXIT_FAILURE);
	}

	if (this->hasExciton(e))
	{
		if (currEn == 1)
		{
			ex1 = nullptr;
			return true;
		}
		if (currEn == 2)
		{
			ex2 = nullptr;
			return true;
		}
	}
	return false; //no exciton present to remove
}

// Checks to see if the exciton that is passes is the exact exciton that already exists in the location.
bool segment::hasExactExciton(shared_ptr<exciton> e)
{
	//If adding the same exciton to the same location, that is self scattering and allowed
	if (e == ex1 || e == ex2)
	{
		return true;
	}
	return false;
}

// calculates the distance between current segment and another segment.
double segment::get_distance(segment &f_seg)
{
	vector<double> diff = vector_sub(point_m , f_seg.point_m);
	double norm = vector_norm(diff);
	return norm;
	
}

// calculates the angle between current segment and another segment.
double segment::get_angle(segment &f_seg)
{
	vector<double> v1 = vector_sub(point1, point2);
	vector<double> v2 = vector_sub(f_seg.point1, f_seg.point2);
	double val = acos(dot_product(v1, v2) / (vector_norm(v1)*vector_norm(v2))); //range 0 to pi
	if (val <= M_PI / 2.0)
	{
		return val;
	}
	return (M_PI - val);
}
