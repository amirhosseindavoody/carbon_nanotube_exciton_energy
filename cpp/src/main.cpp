#include <iostream>
#include <ctime>
#include <armadillo>

#include "cnt.h"

int main(int argc, char *argv[])
{

	std::clock_t start = std::clock();

	std::time_t start_time = std::time(nullptr);
	std::cout << "\nstart time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	cnt m_cnt;

	m_cnt.process_command_line_args(argc, argv);
	m_cnt.geometry();
	m_cnt.electron_full();

	// // cnts[0].dielectric();
	// cnts[0].coulomb_int();

	std::clock_t end = std::clock();

	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;

	return 0;
}
