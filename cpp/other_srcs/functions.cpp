#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include "dirent.h"
#include <list>
#include "CNT.h"
#include <vector>
#include <memory>
#include "tableElem.h"
#include <locale>
#include <math.h>
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include <regex>
#include <omp.h>
#include <stdint.h>

#include "functions.h"

using namespace math_functions;



// Removes excitons from the list if they are in the exit contact and inject excitons into the inContact if there are not enough excitons
void updateExcitonList(int numExcitonsAtCont, vector<exciton> &excitons, vector<int> &currCount, vector<shared_ptr<segment>> inContact)
{
	// Adding excitons to the injection contact
	int numExAdd = numExcitonsAtCont - currCount[0]; //The number of excitons to be added to the injection contact
	if (numExAdd > 0)
	{
		for (int i = 0; i<numExAdd; i++)
		{
			exciton curr_exciton;
			excitons.push_back(curr_exciton); //initialize exciton at end of exciton list
			injectExciton(curr_exciton, inContact); //last element is the exciton to inject
		}
	}

	//Removing excitons from the exit contact
	int numExRem = currCount.back(); //Last element in the count list = num of excitons to remove
	if (numExRem > 0)
	{
		for (int i = 0; i < excitons.size(); i++)
		{
			if (excitons[i].isAtOutContact())
			{
				excitons.erase(excitons.begin() + i);
				i--; // this is to make sure that the new exciton object with now index equal to "i" is not at the out_contact.
			}
		}
	}
}


// Writes the pertinent information of the current state of simulation to file
void writeStateToFile(ofstream &file, vector<int> &currCount, double time)
{
	//Time followed by the counts
	file << std::scientific << time << "     "; 
	for (uint32_t i = 1; i < currCount.size(); i++)
	{
		file << "     " << std::scientific << currCount[i];
	}
	file << endl;
}



// Takes the current exciton and adds it to the count of excitons in a certain region. This is called once per exciton per time step.
void markCurrentExcitonPosition(vector<CNT> &cnt_list, exciton &curr_exciton, vector<int> &currCount, vector<double> &regionBdr)
{
	//add count to the currCount vector in the location corresponding to the exciton region location
	int cnt_idx = curr_exciton.getCNTidx();
	int seg_idx = curr_exciton.getSegidx();
	CNT &curr_cnt = cnt_list[cnt_idx];
	segment &curr_seg = curr_cnt.segments[seg_idx];

	int i = get_index(regionBdr, curr_seg.point_m[0]);
	currCount[i]++;
}


// Checks to see if an exciton has moved into the out contact
bool hasMovedToOutContact(exciton &curr_exciton, vector<double> &regionBdr, vector<CNT> &cnt_list)
{
	//If x component is greater than the second to last index of regionBdr, then it is in the output
	int cnt_idx = curr_exciton.getCNTidx();
	int seg_idx = curr_exciton.getSegidx();
	CNT &curr_cnt = cnt_list[cnt_idx];
	segment &curr_seg = curr_cnt.segments[seg_idx];

	if ((curr_seg.point_m[0]) >= (regionBdr[regionBdr.size() - 2]))
	{
		return true;
	}
	else
	{
		return false;
	}
}



// Places the specified exciton into the input contact
void injectExciton(exciton &curr_exciton, vector<shared_ptr<segment>> &in_contact)
{
	//randomly set the energy of the exciton
	int energy = static_cast<int>(round(getRand(false))+1);
	curr_exciton.setEnergy(energy);

	//choose a destination segment
	int contact_idx = static_cast<int>(rand() % in_contact.size());
	segment &injected_seg = *(in_contact[contact_idx]);

	//The self scattering table will have the correct indices for the current segment
	vector<tableElem> &tbl = injected_seg.tbl;

	//Set the exciton indices to the current segment
	int cnt_idx = (tbl[tbl.size()-1]).getTubeidx();
	int seg_idx = (tbl[tbl.size()-1]).getSegidx();
	curr_exciton.setCNTidx(cnt_idx);
	curr_exciton.setSegidx(seg_idx);
	curr_exciton.setAtOutContact(false); //initialized out contact boolean
}



// Assigns the specified exciton to the next state in the simulation.
void assignNextState(vector<CNT> &cnt_list, exciton &curr_exciton, double gamma, vector<double> &regionBdr)
{

	int cnt_idx = curr_exciton.getCNTidx();
	int seg_idx = curr_exciton.getSegidx();

	CNT &curr_cnt = cnt_list[cnt_idx];
	segment &curr_segment = curr_cnt.segments[seg_idx];
	
	int tblIdx = get_index(curr_segment.rateVec, getRand(false)*gamma);
	tableElem &tbl = curr_segment.tbl[tblIdx];
	
	curr_exciton.setCNTidx(tbl.getTubeidx());
	curr_exciton.setSegidx(tbl.getSegidx());

	bool moved_to_out_contact = hasMovedToOutContact(curr_exciton, regionBdr, cnt_list);
	curr_exciton.setAtOutContact(moved_to_out_contact);
}


// Adds to each segments' rate vector the self scattering component of the simulation
void add_self_scattering(vector<shared_ptr<segment>> &seg_list, double max_rate)
{
	for (int i=0; i<seg_list.size(); i++)
	{
		segment &seg = *(seg_list[i]);
		double self_rate = max_rate - seg.rateVec.back();
		if(self_rate>0)
		{
			int cnt_idx = seg.cnt_idx;
			int seg_idx = seg.seg_idx;
			seg.tbl.push_back(tableElem(1.0, 0.0, self_rate, cnt_idx, seg_idx))
			seg.rateVec.push_back(max_rate);
		}
	}
}



// Gets a random number between 0 and 1
double getRand(bool excludeZero)
{
	if (excludeZero)
	{
		int r;
		while ((r = rand()) == 0){}
		return static_cast<double>(r) / static_cast<double>(RAND_MAX);
	}
	return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}


// calculates the table of scattering rates from a segment to all possible final segments.
// it also returns the total scattering rates out of that element.
double make_rate_table(vector<shared_ptr<segment>> &seg_list, segment &seg, double max_dist)
{
	double rate = 0;

	for (int i = 0; i<seg_list.size(); i++)
	{
		segment &f_segment = *(seg_list[i]);

		double r = seg.get_distance(f_segment);
		double theta = seg.get_angle(f_segment);

		//Check if within range
		if ((r > 0)&&(r <= max_dist)) /////// Building TABLE /////
		{
			double g = 6.4000e+19; //First draft estimate
			seg.tbl.push_back(tableElem(r, theta, g, seg.cnt_idx, seg.seg_idx)); //tbl initialized in CNT::calculateSegments
			rate += (seg.tbl.back()).getRate();
			seg.rateVec.push_back(rate);//tbl initialized in CNT::calculateSegments
		}	
	}

	return rate;
}


// Converts numbers with some units to angstroms
double convert_units(string unit, double val)
{
	if (unit.compare("mm") == 0 || unit.compare("millimeter") == 0)
	{
		return val * 10000000;
	}
	else if (unit.compare("um") == 0 || unit.compare("micrometer") == 0)
	{
		return val * 10000;
	}
	else if (unit.compare("nm") == 0 || unit.compare("nanometer") == 0)
	{
		return val * 10;
	}
	else if (unit.compare("pm") == 0 || unit.compare("picometer") == 0)
	{
		return val * .01;
	}
	else if (unit.compare("A") == 0 || unit.compare("angstrom") == 0)
	{
		return val;
	}
	else
	{
		cout << "error in converting units!!!" << endl;
		exit(EXIT_FAILURE);
	}
}

// Initialized the random number generator by providing a seed value from the time.
void init_random_number_generator()
{
	time_t seconds;
	time(&seconds); //assign time from clock
	srand(static_cast<int>(seconds));
}