#include <iostream>
#include <ctime>

#include "file_management.h"
#include "cnt.h"
#include "constants.h"
#include "write_log.h"
#include "nr3.h"

int main(int argc, char *argv[])
{

	std::clock_t start = std::clock();

	std::time_t start_time = std::time(nullptr);
	std::cout << std::endl << "start time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	vector<cnt> cnts;


	if (argc != 2)
	{
		std::cout << "input directory must be entered as an argument!!!" << std::endl;;
		std::exit(EXIT_FAILURE);
	}

	file_management fm;
	fm.set_input_directory(argv[1]);
	fm.parse_xml(cnts);
	fm.change_working_directory(fm.get_output_directory().c_str());


	cnts[0].geometry();
	cnts[0].electron();

	double a= 2.0;
	double b = nr::sqr(a);

	std::cout << "a = " << a << " , b = " << b << "\n";

	// int *ip = new int[10];
	// int &ir = *ip;

	// for (int i=0; i<10; i++) std::cout << ip[i] << "  ";
	// std::cout << std::endl;
	// for (int i=0; i<10; i++) std::cout << (&ir)[i] << "  ";
	// std::cout << std::endl;


	// // create cnt objects and load their information
	// vector<CNT> cnt_list;
	// for (int i=0; i < file_manager.file_list.size(); i++)
	// {
	// 	string file_name = file_manager.file_list[i];
	// 	string path = file_manager.get_input_directory();
	// 	double segment_length = sim.segment_length;

	// 	cnt_list.push_back(CNT(file_name, path, segment_length));
	// }

	// // create a list of segment pointers that has a list of all segments
	// vector<shared_ptr<segment>> seg_list(0);
	// for (int i=0; i<cnt_list.size(); i++)
	// {
	// 	for (int j=0; j<cnt_list[i].segments.size(); j++)
	// 	{
	// 		cnt_list[i].segments[j].cnt_idx = i;
	// 		cnt_list[i].segments[j].seg_idx = j;
	// 		shared_ptr<segment> seg_ptr = make_shared<segment>(cnt_list[i].segments[j]);
	// 		seg_list.push_back(seg_ptr);
	// 	}
	// }
	
	// {
	// 	string log_input = "number of segments = " + to_string(seg_list.size());
	// 	write_log(log_input);
	// }

	// sim.rmax = sqrt(pow(sim.rmax, 2) + pow(ymax,2));
	// // ********************************************************************************************
	// // get some statistics about cnt network
	// // ********************************************************************************************
	// segment_x_distribution(cnt_list, sim);
	// segment_angle_distance_distribution(seg_list, sim);

	// // ********************************************************************************************
	// // create list of contact segments
	// // ********************************************************************************************
	// int num_regions = static_cast<int>(sim.xdim / sim.region_length_min); //number of regions in the simulation
	// vector<double> region_boundaries = linspace(-sim.xdim/2.0, sim.xdim/2, num_regions); //The boundary of the rgions in the x direction
	// vector<shared_ptr<segment>> in_contact(0); //List of segments in the first region, which is used as a input contact
	// vector<shared_ptr<segment>> out_contact(0); //List of segments in the first region, which is used as a input contact

	// for (int i=0; i<seg_list.size(); i++)
	// {
	// 	segment &curr_segment = *(seg_list[i]);

	// 	int region_idx = get_index(region_boundaries, curr_segment.point_m[0]);

	// 	if (region_idx == 0)
	// 	{
	// 		in_contact.push_back(make_shared<segment>(curr_segment));
	// 	}
	// 	else if((region_idx == region_boundaries.size()-2) && (region_idx == region_boundaries.size()-1))
	// 	{
	// 		out_contact.push_back(make_shared<segment>(curr_segment));
	// 	}
	// }

	// // ********************************************************************************************
	// // build scattering rate tables for each segment
	// // ********************************************************************************************
	// double max_rate = 0; //maximum total scattering rate from all segments.
	// for (int i=0; i<seg_list.size(); i++)
	// {
	// 	segment &curr_segment = *(seg_list[i]);
	// 	double total_rate = make_rate_table(seg_list, curr_segment, sim.maximum_distance);
	// 	if (total_rate > max_rate)
	// 	{
	// 		max_rate = total_rate;
	// 	}		
	// }

	// // write the maximum scattering rate in the simulation
	// {
	// 	stringstream log_input;
	// 	log_input << "maximum scattering rate = " << std::scientific << max_rate;
	// 	write_log(log_input.str());
	// }

	// // add self-scattering rates to the scattering table
	// add_self_scattering(seg_list, max_rate);

	
	// // ********************************************************************************************
	// // populate segments with excitons
	// // ********************************************************************************************
	// vector<exciton> excitons(sim.number_of_excitons);
	// for (int i = 0; i<excitons.size(); i++)
	// {
	// 	injectExciton(excitons[i], in_contact);
	// }


	// //////////////////////////////////// TIME STEPS ///////////////////////////////////
	// double deltaT = (1 / max_rate)*sim.tfac; //time steps at which statistics are calculated
	// double Tmax = deltaT * sim.number_of_steps; //maximum simulation time
	// double T = 0; //Current simulation time, also time at which next stats will be calculated

	// vector<int> currCount(num_regions); //The count of excitons in each region of the CNT mesh

	// //File output initializations
	// string excitonDistFileName = file_manager.output_directory + "excitonDist.csv";
	// ofstream excitonDistFile;
	// excitonDistFile.open(excitonDistFileName);

	// //For automatic simulation end there will be a comparison of differences to some sim.threshold
	// // If the differences, divided by the maximum current differences are < sim.threshold for more
	// // the sim.num_to_finish then the simulation will end.
	// int numInARow = 0; //number of quotients that are below sim.threshold in a row
	// double difference = 0;//difference between average of sim.num_to_check time step average num exciton points
	// double maxDiff = 0; //The maximum difference between points used for difference.
	// double prevAve =  sim.number_of_excitons; //average for last sim.num_to_check time steps. Starts at init num of excitons
	// double currAve = 0; //Average for current sim.num_to_check time steps
	// int timeSteps = 0; //current count of number of time steps
	// bool simDone = false; //boolean symbolizing a finished simulation

	// // ********************************************************************************************
	// // monte carlo stepping loop
	// // ********************************************************************************************
	// while (T <= Tmax && !simDone) //iterates over time
	// {
	// 	T += deltaT;

	// 	//reset the exciton count after each time step
	// 	vector_clear(currCount, 0);

	// 	for (int exNum = 0; exNum < excitons.size(); exNum++)
	// 	{
	// 		exciton &curr_exciton = excitons[exNum];
	// 		/*There is a change that the previous tr that was calculated was so long that it
	// 		not only has extra time in the next deltaT, but it skips it completely. In this case
	// 		the exciton movement is skipped all together and its extra time is decreased by deltaT
	// 		until the exciton can move again. Otherwise the code runs as usual.
	// 		*/
	// 		double extraT = curr_exciton.getTExtra();
	// 		if (extraT > deltaT)
	// 		{
	// 			curr_exciton.setTExtra(extraT - deltaT);
	// 			markCurrentExcitonPosition(cnt_list, curr_exciton, currCount, region_boundaries);
	// 		}
	// 		else
	// 		{
	// 			double tr_tot = extraT; //the sum of all tr's in the current deltaT time step
	// 			/*give time to excitons that have no extra time. This happens only when excitons
	// 			are just injected.*/
	// 			if (extraT == 0)
	// 			{
	// 				tr_tot += -(1 / max_rate)*log(getRand(true));
	// 			}
	// 			while (tr_tot <= deltaT)
	// 			{
	// 				//choose new state
	// 				assignNextState(cnt_list, curr_exciton, max_rate, region_boundaries);
	// 				tr_tot += -(1 / max_rate)*log(getRand(true)); // add the tr calculation to current time for individual particle
	// 			}
	// 			//Recording distribution
	// 			markCurrentExcitonPosition(cnt_list, curr_exciton, currCount, region_boundaries);
	// 			//record the time past deltaT that the assign next state will cover
	// 			curr_exciton.setTExtra(tr_tot - deltaT);
	// 		}
	// 	}

	// 	if (sim.auto_complete)
	// 	{
	// 		timeSteps++;
	// 		currAve += excitons.size();
	// 		if (sim.num_to_check == timeSteps)
	// 		{
	// 			currAve /= sim.num_to_check; //calculate average
				
	// 			//check for max difference
	// 			if ((difference = currAve - prevAve) > maxDiff)
	// 			{
	// 				maxDiff = difference;
	// 			}
				
	// 			if (difference / maxDiff < sim.threshold)
	// 			{
	// 				if (++numInARow >= sim.num_to_finish)
	// 				{
	// 					simDone = true;
	// 				}
	// 			}
	// 			else
	// 			{
	// 				numInARow = 0;
	// 			}
	// 			prevAve = currAve;
	// 			timeSteps = 0; //reset the time counter
	// 		}
	// 	}

	// 	//Output count vector to file since we want results after each time step.
	// 	writeStateToFile(excitonDistFile, currCount, T);

	// 	//Update Exciton List for injection and exit contact
	// 	updateExcitonList(sim.number_of_excitons, excitons, currCount, in_contact);
		
	// }

	// //Close files finish program
	// excitonDistFile.close();

	return 0;
}