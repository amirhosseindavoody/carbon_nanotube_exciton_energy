#include <iostream>
#include <ctime>
#include <armadillo>

#include "cnt.h"
#include "exciton_transfer.h"

int main(int argc, char *argv[])
{
	std::clock_t start = std::clock();

	std::time_t start_time = std::time(nullptr);
	std::cout << "\nstart time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	cnt m_cnt;

	m_cnt.process_command_line_args(argc, argv);
	m_cnt.calculate_exciton_dispersion();

	exciton_transfer ex_transfer(m_cnt, m_cnt);
	ex_transfer.save_Q_matrix_element(0,0);
	ex_transfer.save_J_matrix_element(0,0);
	// ex_transfer.first_order();

	std::clock_t end = std::clock();

	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;

	return 0;
}
