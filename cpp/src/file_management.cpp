#include <stdio.h>
#include <iostream>
#include <string>
#include <dirent.h>
#include <regex>
#include <vector>
#include <sys/stat.h> // used for checking status of folders and files and see if they exist
#include <unistd.h> // used for changing the working directory in the program: chdir

#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"

#include "file_management.h"
#include "cnt.h"

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
		std::cout << input_directory << std::endl;
		std::cout << "error: incorrect result directory path!!!" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	// check if there is '/' at the end of the input directory
	char last_char = input_directory.at(input_directory.size()-1);
	if (last_char != '/')
	{
		input_directory.push_back('/');
	}
	std::cout << "input directory set to: " << input_directory << std::endl;

}

// this function sets the input directory and checks if it exists and accessible.
void file_management::set_output_directory(const char *output)
{
	output_directory = output;

	// check if there is '/' at the end of the input directory
	char last_char = output_directory.at(output_directory.size()-1);
	if (last_char != '/')
	{
		output_directory.push_back('/');
	}

	// // check if folder or file exists.
	// struct stat buf;
	// if (stat(output_directory.c_str(),&buf) != 0)
	// {
	// 	std::cout << "error: incorrect output directory path!!!" << std::endl;;
	// 	std::exit(EXIT_FAILURE);
	// }

	std::cout << "output directory set to: " << output_directory << std::endl;
}

// this function gives back the input directory
std::string file_management::get_input_directory()
{
	return input_directory;
}

// this function gives back the output directory
std::string file_management::get_output_directory()
{
	return output_directory;
}

// parse the input xml file and creat the cnt objects with input properties
void file_management::parse_xml(std::vector<cnt> &cnts)
{
	std::string filename = get_input_directory() + "input.xml";

	rapidxml::xml_document<> doc; //create xml object
	rapidxml::file<> xmlFile(filename.c_str()); //open file
	doc.parse<0>(xmlFile.data()); //parse contents of file

	rapidxml::xml_node<>* node = doc.first_node(); //gets the node "Document" or the root node

	node = node->first_node("output_directory");
	std::string outdir = node->value();
	set_output_directory(outdir.c_str());

	node = node->next_sibling("cnt");
	std::string name = node->first_node("name")->value();
	int n = atoi(node->first_node("chirality")->first_node("n")->value());
	int m = atoi(node->first_node("chirality")->first_node("m")->value());
	int length = atoi(node->first_node("length")->first_node("value")->value());

	// std::cout << "outdir = " << outdir << " , name = " << name << " , n = " << n << " , m = " << m << " , length = " << length << std::endl;
	cnts.push_back(cnt(name, n, m, length));

	return;
}


void file_management::change_working_directory(const char *path)
{
	// check if folder or file exists.
	struct stat buf;
	if (stat(path,&buf) == 0)
	{
		std::stringstream command;
		command << "rm -f -r " << path;
		system(command.str().c_str());
	}

	std::stringstream command;
	command << "mkdir " << path;
	system(command.str().c_str());

	int rc = chdir(path);
	if (rc < 0)
	{
		std::cout << "unable to change the output directory!!!" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	const int max_path_length = 256;
	char buffer[max_path_length];
	char *current_path = getcwd(buffer, max_path_length);
	std::cout << "current working directory: "<< current_path << std::endl;
}