#pragma once

#include <Eigen>
#include <memory>

#include "tableElem.h"
#include "exciton.h"
#include "math_functions.h"

using namespace std;
using namespace Eigen;

using namespace math_functions;
class tableElem;

//Stores the position information about 
struct segment
{
	//instance variables
	int segment_number;
	vector<tableElem> tbl;
	vector<double> rateVec; // Vector to store all rates

	vector<double> point1;
	vector<double> point2;
	vector<double> point_m;

	int cnt_idx;
	int seg_idx;
	
	shared_ptr<exciton> ex1; //first energy level for exciton
	shared_ptr<exciton> ex2; //second energy level for exciton

	segment(int number, vector<double> first_point, vector<double> second_point);
	bool hasExciton(shared_ptr<exciton> e);
	bool setExciton(shared_ptr<exciton> e);
	bool removeExciton(shared_ptr<exciton> e);
	bool hasExactExciton(shared_ptr<exciton> e);

	double get_distance(segment &f_seg);
	double get_angle(segment &f_seg);
};
