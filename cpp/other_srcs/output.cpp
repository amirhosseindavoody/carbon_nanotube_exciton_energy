#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <vector>
#include <memory>

#include "simulation_parameters.h"
#include "CNT.h"
#include "output.h"
#include "functions.h"
#include "segment.h"

using namespace std;

// this function calculates and saves the number of cnt segments along the x-axis in the simulation box.
void output::segment_x_distribution(vector<CNT> &cnt_list, simulation_parameters sim)
{
	int num_regions = static_cast<int>(sim.xdim / sim.region_length_min); //number of regions in the simulation
	vector<double> region_boundaries = linspace(-sim.xdim/2.0, sim.xdim/2, num_regions); //The boundary of the rgions in the x direction
	vector<int> seg_per_region(num_regions); //To get dist stats

	// count the number of segments per region
	for (int i=0; i<cnt_list.size(); i++)
	{
		CNT &curr_cnt = cnt_list[i];
		
		for (int j=0; j<curr_cnt.segments.size(); j++)
		{
			segment &curr_segment = curr_cnt.segments[j];

			int regIdx = get_index(region_boundaries, curr_segment.point_m[0]);
			seg_per_region[regIdx]++; //increment the count based on where segment is
		}		
	}

	// save the results
	string file_name = "segment_per_region.dat";
	ofstream out_file;
	out_file.open(file_name);
	for (int i = 0; i < seg_per_region.size(); i++)
	{
		out_file << std::scientific << seg_per_region[i] << endl;
	}
	out_file.close();

}


void output::segment_angle_distance_distribution(vector<shared_ptr<segment>> &seg_list, simulation_parameters sim)
{
	int n_r = 100;
	int n_theta = 100;
	vector<double> r_vec = linspace(0, sim.rmax, n_r);
	vector<double> theta_vec = linspace(0, 90, n_theta);
	vector<vector<double>> heat_map(n_r, vector<double>(n_theta));
	
	//loop over CNTs
	for (int i=0; i<seg_list.size(); i++)
	{
		segment &seg1 = *(seg_list[i]);
		
		for (int j=0; j<seg_list.size(); j++)
		{
			segment &seg2 = *(seg_list[j]);

			double r = seg1.get_distance(seg2);
			double theta = seg1.get_angle(seg2);

			int i_r = get_index(r_vec, r);
			int i_theta = get_index(theta_vec, theta);

			heat_map[i_r][i_theta]++;
		}		
	}

	// save the generated heat_map
	string file_name = "heat_map.dat";
	ofstream out_file;
	out_file.open(file_name);
	for (int i = 0; i < r_vec.size(); i++)
	{
		for (int j = 0; j < theta_vec.size(); j++)
		{
			out_file << std::scientific << heat_map[i][j] << "         ";
		}
		out_file << endl;
	}
	out_file.close();

}