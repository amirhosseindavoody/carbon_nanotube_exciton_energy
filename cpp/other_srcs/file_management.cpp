#include <stdio.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <regex>
#include <list>
#include <vector>
#include <memory>
#include <unistd.h> // used for changing the working directory in the program: chdir

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "file_management.h"
#include "simulation_parameters.h"
#include "functions.h"

using namespace std;

file_management:: file_management()
{
}

file_management:: ~file_management()
{
}

// this function sets the input directory and checks if it exists and accessible.
void file_management::set_input_directory(const char *input)
{
	input_directory = input;

	// check if folder or file exists.
	struct stat buf;
	if (stat(input_directory.c_str(),&buf) != 0)
	{
		cout << "Error: incorrect result directory path!!!" << endl;;
		exit(EXIT_FAILURE);
	}

	// check if there is '/' at the end of the input directory
	char last_char = input_directory.at(input_directory.size()-1);
	if (last_char != '/')
	{
		input_directory.push_back('/');
	}
	cout << "input directory: " << input_directory << endl;

}

// this function gives back the input directory
string file_management::get_input_directory()
{
	return input_directory;
}

// This function creates a list of files that contain information of cnt geometry and location from the bullet physics simulation.
void file_management::build_cnt_file_list()
{
	DIR *resDir;
	struct dirent *ent = nullptr;

	regex rgx("CNT_Num_\\d+\\.csv"); //files we want to look through

	//Check if folder can be opened - should work due to above checks
	if ((resDir = opendir(input_directory.c_str())) != nullptr)
	{
		//iterate over all of the real files
		while ((ent = readdir(resDir)) != nullptr)
		{
			smatch matches; //match_results for string objects
			string tmps =  string(ent->d_name);
			regex_search(tmps, matches, rgx);
			if (!matches.empty())
			{
				file_list.push_back(ent->d_name);
			}
		}
		closedir(resDir); //deletes pointer
	}
	else
	{
		cout << "Could not open input directory!!!" << endl;
		exit(EXIT_FAILURE);
	}
	delete ent;
	ent = nullptr;
}

simulation_parameters file_management::parse_xml()
{
	simulation_parameters sim;

	string inputXMLPath = get_input_directory() + "input.xml";

	rapidxml::xml_document<> doc; //create xml object
	rapidxml::file<> xmlFile(inputXMLPath.c_str()); //open file
	doc.parse<0>(xmlFile.data()); //parse contents of file
	rapidxml::xml_node<>* currNode = doc.first_node(); //gets the node "Document" or the root node


	currNode = currNode->first_node("outputDirectory");
	output_directory = currNode->value();

	cout << "output_directory: " << output_directory << endl;

	// DEVICE DIMENSIONS NODE //
	currNode = currNode->next_sibling("DeviceDimensions"); //Output folder

	sim.xdim = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("xdim")->value()));
	sim.ydim = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("ydim")->value()));
	sim.zdim = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("zdim")->value()));

	if (sim.xdim <= 0 || sim.ydim <= 0 || sim.zdim <= 0)
	{
		cout << "error: must enter positive device dimensions!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END DEVICE DIMENSIONS NODE //

	//x and z dim are the ground dimensions and y is height. Need to make sure bottom corner
	// to top other corner is included in r range. ymax is not defined so this is intermediate
	// value. It is changed after CNT initialization.
	sim.rmax = sqrt(pow(sim.xdim,2)+pow(sim.zdim,2));

	// NUMBER OF EXCITONS NODE //
	currNode = currNode->next_sibling("numberExcitons"); //to number of excitons
	sim.number_of_excitons = atoi(currNode->value());
	if (sim.number_of_excitons <= 0)
	{
		cout << "Error: Must have positive number of excitons." << endl;
		exit(EXIT_FAILURE);
	}
	// END NUMBER OF EXCITONS NODE //


	// REGION LENGTH NODE //
	currNode = currNode->next_sibling("regionLength");
	sim.region_length_min = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
	if (sim.region_length_min <= 0 || sim.region_length_min > sim.xdim)
	{
		cout << "Error: Region length must be positive number!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END REGION LENGTH NODE //

	// SEGMENT LENGTH NODE //
	currNode = currNode->next_sibling("segmentLength");
	sim.segment_length = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
	if (sim.segment_length <= 0)
	{
		cout << "Error: Segment length must be positive number!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END SEGMENT LENGTH NODE //

	// SEGMENT SEPARATION //
	currNode = currNode->next_sibling("segmentSeparation");
	sim.maximum_distance = convert_units(string(currNode->first_node("units")->value()), atof(currNode->first_node("min")->value()));
	if (sim.maximum_distance <= 0)
	{
		cout << "Error: Segment separation must be a positive number!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END SEGMENT SEPARATION //

	// NUMBER OF TIME STEPS //
	currNode = currNode->next_sibling("numberTimeSteps");
	sim.number_of_steps = atoi(currNode->value());
	if (sim.number_of_steps <= 0)
	{
		cout << "Error: Must have positive number time steps.!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END NUMBER OF TIME STEPS //


	// PERCENT FREE FLIGHT TIMES ABOVE DELTA T //
	currNode = currNode->next_sibling("percentFreeFlightTimesAboveDeltaT");
	sim.tfac = atoi(currNode->value()) / 100.0;
	if (sim.tfac <= 1 && sim.tfac > 0)
	{
		sim.tfac = abs(log(sim.tfac));
	}
	else
	{
		cout << "Error: Percent of free flight times above delta T must be greater than 0 and less than or equal to 100!!!" << endl;
		exit(EXIT_FAILURE);
	}
	// END PERCENT FREE FLIGHT TIMES ABOVE DELTA T  //

	// AUTO COMPLETE //
	currNode = currNode->next_sibling("autoComplete")->first_node("enabled");

	if (string(currNode->value()).compare("true") == 0)
	{
		sim.auto_complete = true;
		sim.number_of_steps = 10000000; // set the number of time steps to a very large number so that it would not be reached before the auto complete is reached!

		currNode = currNode->next_sibling("threshold");
		sim.threshold = atof(currNode->value());
		if (sim.threshold < 0.0 || sim.threshold > 1.0)
		{
			cout << "Configuration Error: threshold must be value greater than 0 and less than 1!!!" << endl;
			exit(EXIT_FAILURE);
		}

		currNode = currNode->next_sibling("numBelowThreshold");
		sim.num_to_finish = atoi(currNode->value());
		if (sim.num_to_finish <= 0) //input validation
		{
			cout << "Configuration Error: Number below threshold must be greater than 0!!!" << endl;
			exit(EXIT_FAILURE);
		}

		currNode = currNode->next_sibling("numToAverage");
		sim.num_to_check = atoi(currNode->value());
		if (sim.num_to_check <= 0)
		{
			cout << "Configuration Error: Number to average must be greater than 0!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}
	// END AUTO COMPLETE // 

	return sim;
}


void file_management::change_working_directory(const char *path)
{
	// check if folder or file exists.
	struct stat buf;
	if (stat(path,&buf) == 0)
	{
		stringstream command;
		command << "rm -f -r " << path;
		system(command.str().c_str());
	}

	stringstream command;
	command << "mkdir " << path;
	system(command.str().c_str());

	int rc = chdir(path);
	if (rc < 0)
	{
		cout << "unable to change the output directory!!!" << endl;
		exit(EXIT_FAILURE);
	}

	// const int max_path_length = 256;
	// char buffer[max_path_length];
	// char *current_path = getcwd(buffer, max_path_length);
	// cout << "current working directory: "<< current_path << endl;
}