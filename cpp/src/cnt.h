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

	int _Nu; // number of graphene unit cells in cnt unit cell.

	arma::vec _K1; // cnt reciprocal lattice vector in the circumferencial direction
	arma::vec _K2; // cnt reciprocal lattice vector along the cnt axis
	arma::vec _K2_normed; // normalised _K2 vector (cnt reciprocal lattice vector along the cnt axis) for generating k vector
	arma::vec _dk_l; // delta_k in the longitudinal direction with respect to cnt axis

	arma::mat _pos_a, _pos_b; // position of atoms in A and B sites
	arma::mat _pos_2d, _pos_3d; // position of all atoms in cnt unit cell in 2d and in 3d space
	arma::mat _pos_aa, _pos_ab, _pos_ba, _pos_bb; // distance between atoms and A and B sites which conserves lattice symmetry.

	arma::mat _el_energy_full; // energy of electronic states calculated using the full unit cell (2*Nu atoms)
	arma::cx_cube _el_psi_full; // electronic wave functions corresponding to electronic states using the full unit cell (2*Nu atoms)

	arma::cube _el_energy_redu; // energy of electronic states calculated using the reduced graphene unit cell (2 atoms)
	arma::field<arma::cx_cube> _el_psi_redu; // electronic wave functions corresponding to electronic states using the reduced graphene unit cell (2 atoms)

	arma::cx_vec _epsilon; // static dielectric function


public:
	//constructor
	cnt(){};
	// set the output directory and the output file name
	void process_command_line_args(int argc, char* argv[]);
	// calculates position of atoms and reciprocal lattice vectors
	void geometry();
	// calculate electron dispersion energies using full unit cell (2*Nu atoms)
	void electron_full();
	// calculate electron dispersion energies using the reduced graphene unit cell (2 atoms)
	void electron_reduced();
	// given index i_mu between [0,Nu-1] return the value of mu between [1-Nu/2, Nu/2]
	const double mu(const int& i_mu) const
	{
		return double(i_mu+1-_Nu/2);
	};
	// given index ik shift it by -_nk/2;
	const double k(const int& ik) const
	{
		return double(ik-_nk/2);
	};
	// void dielectric(); // calculate static dielectric function
	// void coulomb_int(); // calculate coulomb interaction matrix elements
};

#endif // end _cnt_h_
