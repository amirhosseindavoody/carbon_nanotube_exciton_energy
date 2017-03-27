
#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>

#include "simulation_parameters.h"

using namespace std;

simulation_parameters::simulation_parameters()
{
	segment_length = 100.0;
	region_length_min = segment_length;
	number_of_steps = 10000;
	maximum_distance = 300.0;
	number_of_excitons = 10;
	auto_complete = false;
	threshold = 0.01;
	num_to_finish = 5;
	num_to_check = 1000;
	tfac = log(0.3);

	rmax = 0;
	xdim = 0;
	ydim = 0;
	zdim = 0;
}