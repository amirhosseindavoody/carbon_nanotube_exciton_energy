#ifndef output_h
#define output_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <vector>
#include <memory>

#include "CNT.h"
#include "simulation_parameters.h"

using namespace std;

namespace output
{
	void segment_x_distribution(vector<CNT> &cnt_list, simulation_parameters sim);
	void segment_angle_distance_distribution(vector<shared_ptr<segment>> &seg_list, simulation_parameters sim);
};


#endif // output_h