#include <iostream>
#include <ctime>

#include "file_management.h"
#include "cnt.h"
#include "constants.h"
#include "write_log.h"
#include "nr.h"

int main(int argc, char *argv[])
{

	std::clock_t start = std::clock();

	std::time_t start_time = std::time(nullptr);
	std::cout << std::endl << "start time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	std::vector<cnt> cnts;


	if (argc != 2)
	{
		std::cout << "input directory must be entered as an argument!!!" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	file_management fm;
	fm.set_input_directory(argv[1]);
	fm.parse_xml(cnts);
	fm.change_working_directory(fm.get_output_directory().c_str());


	cnts[0].geometry();
	cnts[0].electron();
	// cnts[0].dielectric();
	cnts[0].coulomb_int();

	std::clock_t end = std::clock();

	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;

	return 0;
}