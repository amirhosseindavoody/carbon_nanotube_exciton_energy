#ifndef _cnt_h_
#define _cnt_h_

#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>

#include "constants.h"
#include "small.h"
#include "../lib/rapidxml.hpp"
#include "../lib/rapidxml_utils.hpp"

class cnt
{
private:

	std::experimental::filesystem::directory_entry _directory; // this is the address of the directory that the cnt data is stored in
	std::string _name; //cnt name

	const double _a_cc = 1.42e-10; // carbon-carbon distance [meters]
	const double _a_l = sqrt(3.0)*_a_cc; // graphene lattice constants [meters]

	const double _e2p = 0.0; // tight binding constants
	const double _t0 = 2.7 * constants::eV; // tight binding constants
	const double _s0 = 0.0; // tight binding constants

	const double _Upp = 11.3 * constants::eV; // constant used in the Ohno potential
	const double _kappa = 2.0; // dielectric constant due to core electrons and surrounding environment

	int _n, _m; // chirailty parameters
	int _number_of_cnt_unit_cells; // length of cnt in units of cnt unit cell.
	int _nk; // number of k vector elements corresponding to the length of the cnt.

	arma::vec _a1, _a2; // real space lattice primitive vectors
	arma::vec _b1, _b2; // reciprocal lattice primitive vectors
	arma::vec _aCC_vec; // vector between two neighboring carbon atoms
	arma::vec _ch_vec; // chirality vector

	double _ch_len; // length of the chirality vector
	double _radius; // radius of the cnt

	arma::vec _t_vec; // translational vector for cnt unit cell in unrolled graphene sheet
	arma::vec _t_vec_3d; // translational vector for cnt unit cell in rolled graphene sheed (3d)
	//
	// int Nu; // number of graphene unit cells in cnt unit cell.
	//
	// nr::vec_doub K1; // cnt reciprocal lattice vector in the circumferencial direction
	// nr::vec_doub K2; // cnt reciprocal lattice vector along the cnt axis
	// nr::vec_doub dk_l; // delta_k in the longitudinal direction with respect to cnt axis
	//
	// nr::mat_doub pos_a, pos_b; // position of atoms in A and B sites
	// nr::mat_doub pos_2d, pos_3d; // position of all atoms in cnt unit cell in 2d and in 3d space
	// // nr::mat_doub pos_aa, pos_ab, pos_ba, pos_bb; // distance between atoms and A and B sites which conserves lattice symmetry.
	//
	// nr::mat_doub el_energy; // energy of electronic states
	// nr::mat3d_complex el_psi; // electronic wave functions corresponding to electronic states
	//
	// nr::vec_complex epsilon; // static dielectric function


public:
	//constructor
	cnt(){};

	// cnt(const std::string &in_name, const int in_n, const int in_m, const int in_length);

	// set the output directory and the output file name
	void process_command_line_args(int argc, char* argv[])
	{
		namespace fs = std::experimental::filesystem;

		// first find the xml input file and open it
		fs::directory_entry xml_file;
		std::cout << "current path is " << fs::current_path() << std::endl;

		if (argc <= 1)
		{
			xml_file.assign("input.xml");
		}
		else
		{
			xml_file.assign(argv[1]);
		}

		if(fs::exists(xml_file))
		{
			std::cout << "input xml file found: " << xml_file.path() << std::endl;
		}
		else
		{
			std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
			std::exit(1);
		}

		if (!fs::is_regular_file(xml_file))
		{
			std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
			std::exit(1);
		}
		std::cout << std::endl;

		rapidxml::file<> xmlFile(xml_file.path().c_str()); //open file
		rapidxml::xml_document<> doc; //create xml object
		doc.parse<0>(xmlFile.data()); //parse contents of file
		rapidxml::xml_node<>* curr_node = doc.first_node(); //gets the node "Document" or the root nodes
		curr_node = curr_node->first_node();

		// get the sibling node with name sibling_name
		auto get_sibling = [](rapidxml::xml_node<>* node, const std::string& sibling_name)
		{
			if (node->name() == sibling_name)
			{
				return node;
			}
			auto next_node = node->next_sibling(sibling_name.c_str());
			if (next_node == 0)
			{
				next_node = node->previous_sibling(sibling_name.c_str());
				if (next_node == 0)
				{
					std::cout << "sibling not found: " << sibling_name.c_str() << std::endl;
					std::exit(1);
				}
			}
			return next_node;
		};

		// get name cnt name
		{
			curr_node = get_sibling(curr_node, "name");
			_name = trim_copy(curr_node->value());
			std::cout << "cnt name: '" << _name << "'\n";
		}

		// set the output_directory
		{
			curr_node = get_sibling(curr_node,"output_directory");

			std::string attr = curr_node->first_attribute("type")->value();
			// std::string path = trim_copy(curr_node->value());
			fs::path path = trim_copy(curr_node->value());
			if (attr == "absolute")
			{
				std::cout << "absolute directory format used!\n";
			}

			_directory.assign(path);
			std::cout << "output_directory: " << _directory.path() << std::endl;

			if (not fs::exists(_directory.path()))
			{
				std::cout << "warning: output directory does NOT exist!!!" << std::endl;
				std::cout << "output directory: " << _directory.path() << std::endl;
				fs::create_directories(_directory.path());
			}

			if (fs::is_directory(_directory.path()))
			{
				if (not fs::is_empty(_directory.path()))
				{
					std::cout << "warning: output directory is NOT empty!!!" << std::endl;
					std::cout << "output directory: " << _directory.path() << std::endl;
					std::cout << "deleting the existing directory!!!" << std::endl;
					fs::remove_all(_directory.path());
					fs::create_directories(_directory.path());
				}
			}
			else
			{
				std::cout << "error: output path is NOT a directory!!!" << std::endl;
				std::cout << "output path: " << _directory.path() << std::endl;
				std::exit(EXIT_FAILURE);
			}
		}

		// read chirality
		{
			curr_node = get_sibling(curr_node,"chirality");
			std::string chirality = curr_node->value();
			_n = std::stoi(chirality.substr(0,chirality.find(",")));
			_m = std::stoi(chirality.substr(chirality.find(",")+1));
			std::cout << "chirality: (" << _n << "," << _m << ")\n";
		}

		// length of cnt
		{
			curr_node = get_sibling(curr_node,"length");
			std::string units = curr_node->first_attribute("units")->value();
			if (units != "cnt_unit_cell")
			{
				std::cout << "cnt length is not in units '" << units << "'\n";
				std::exit(1);
			}
			_number_of_cnt_unit_cells = std::stoi(curr_node->value());
			std::cout << "length of cnt: " << _number_of_cnt_unit_cells << " unit cells.\n";
		}


	};
	// initialize
	void init();

	// calculates position of atoms
	void geometry();
	// void electron(); // electron dispersion energies
	// void dielectric(); // calculate static dielectric function
	// void coulomb_int(); // calculate coulomb interaction matrix elements
};

#endif // end _cnt_h_
