#ifndef file_management_h
#define file_management_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>

#include <boost/property_tree/xml_parser.hpp>

#include "cnt.h"

class file_management
{
	private:
		std::string input_directory; // this is the address of the input folder.
		std::string output_directory; // this is the address of the output folder that the simulation results would be saved into.

	public:
		
		file_management();
		~file_management();
		void set_input_directory(const char *input); // this function sets the input directory and checks if it exists and accessible.
		void set_output_directory(const char *output); // this function sets the output directory and checks if it exists and accessible.
		std::string get_input_directory(); // this function gives back the input directory
		std::string get_output_directory(); // this function gives back the output directory
		void parse_xml(std::vector<cnt> &cnts); // this function parses the input xml file and puts the simulation parameters in the simulation parameter object "sim".
		void change_working_directory(const char *path);		
};


#endif // file_management_h