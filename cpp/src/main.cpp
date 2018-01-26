#include <iostream>
#include <ctime>
#include <armadillo>

#include "cnt.h"
#include "exciton_transfer.h"
#include "constants.h"
#include "../lib/json.hpp"

int main(int argc, char *argv[])
{
	using json = nlohmann::json;

	std::string filename;
	if (argc <= 1){
		filename = "input.json";
	} else {
		filename = argv[1];
	}

	// read a JSON file
	std::ifstream input_file(filename.c_str());
	json j;
	input_file >> j;

	std::string parent_directory = j["cnts"]["directory"];
	j["cnts"].erase("directory");

	std::vector<cnt> cnts;
	for (const auto& j_cnt: j["cnts"])
	{
		cnts.emplace_back(cnt(j_cnt,parent_directory));
	};


	// std::exit(0);

	for (auto& cnt: cnts)
	{
		cnt.calculate_exciton_dispersion();
	}
   

	std::exit(0);




	std::time_t start_time = std::time(nullptr);
	std::cout << "\nstart time:" << std::endl << std::asctime(std::localtime(&start_time)) << std::endl;

	// cnt m_cnt;

	// m_cnt.process_command_line_args(argc, argv);
	// m_cnt.calculate_exciton_dispersion();

	// exciton_transfer ex_transfer(m_cnt, m_cnt);
	// ex_transfer.save_J_matrix_element(0,0);
	// ex_transfer.first_order(1.5e-9, {0.,0.}, 0, true);
	// ex_transfer.calculate_first_order_vs_angle(1.5e-9, {0.,0.});
	// ex_transfer.calculate_first_order_vs_zshift({0.,0.}, 0);
	// ex_transfer.calculate_first_order_vs_axis_shift(1.5e-9, constants::pi/2);


	std::time_t end_time = std::time(nullptr);
	std::cout << std::endl << "end time:" << std::endl << std::asctime(std::localtime(&end_time));
	std::cout << "runtime: " << std::difftime(end_time,start_time) << " seconds" << std::endl << std::endl;

	return 0;
}
