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
	const double _a_l = std::sqrt(float(3.0))*_a_cc; // graphene lattice constants [meters]

	const double _e2p = 0.0; // tight binding constants
	const double _t0 = 2.7 * constants::eV; // tight binding constants
	const double _s0 = 0.0; // tight binding constants

	const double _Upp = 11.3 * constants::eV; // constant used in the Ohno potential
	const double _kappa = 2.0; // dielectric constant due to core electrons and surrounding environment

	int _n, _m; // chirailty parameters
	int _number_of_cnt_unit_cells; // length of cnt in units of cnt unit cell.
	int _nk_K1; // number of k vector elements in the K1-extended representation calculated based on the length of the cnt.

	arma::vec _a1, _a2; // real space lattice primitive vectors
	arma::vec _b1, _b2; // reciprocal lattice primitive vectors
	arma::vec _aCC_vec; // vector between two neighboring carbon atoms
	arma::vec _ch_vec; // chirality vector

	double _ch_len; // length of the chirality vector
	double _radius; // radius of the cnt

	int _t1, _t2; // translation vector elemens in terms of _a1 and _a2
	int _M, _Q; // cnt parameters for K2-extended representation, Q is the number of cutting line in this representation

	arma::vec _t_vec; // translational vector for cnt unit cell in unrolled graphene sheet
	arma::vec _t_vec_3d; // translational vector for cnt unit cell in rolled graphene sheed (3d)

	int _Nu; // number of graphene unit cells in cnt unit cell.

	arma::vec _K1; // cnt reciprocal lattice vector in the circumferencial direction
	arma::vec _K2; // cnt reciprocal lattice vector along the cnt axis
	arma::vec _K2_normed; // normalised _K2 vector (cnt reciprocal lattice vector along the cnt axis) for generating k vector
	arma::vec _dk_l; // delta_k in the longitudinal direction with respect to cnt axis

	arma::mat _pos_a, _pos_b; // position of atoms in A and B sites
	arma::mat _pos_2d, _pos_3d; // position of all atoms in cnt unit cell in 2d and in 3d space

	arma::mat _el_energy_full; // energy of electronic states calculated using the full unit cell (2*Nu atoms)
	arma::cx_cube _el_psi_full; // electronic wave functions corresponding to electronic states using the full unit cell (2*Nu atoms)

	arma::cube _el_energy_redu; // energy of electronic states calculated using the reduced graphene unit cell (2 atoms)
	arma::field<arma::cx_cube> _el_psi_redu; // electronic wave functions corresponding to electronic states using the reduced graphene unit cell (2 atoms)

	int _ik_min_K2, _ik_max_K2; // limits of ik-vector in the K2-extended representation
	int _nk_K2; // number of ik vector elements in the K2-extended representation
	int _mu_min_K2, _mu_max_K2; // limits of mu in the K2-extended representation
	int _n_mu_K2; // number of cutting lines in the K2-extended representation
	arma::cube _el_energy_K2; // energy of electronic states in K2-extended rep. calculated using the reduced graphene unit cell (2 atoms)
	arma::field<arma::cx_cube> _el_psi_K2; // electronic wave functions in K2-extended rep. corresponding to electronic states using the reduced graphene unit cell (2 atoms)
	std::vector<std::array<std::array<unsigned int, 2>, 2>> _valleys_K2; // index of valleys in K2-extended representation
	std::vector<std::vector<std::array<int,2>>> _relev_ik_range;

	// // struct to bundle information of electronic energy state
	// struct el_energy_struct
	// {
	// 	arma::cube e; // energy of electronic states calculated using the reduced graphene unit cell (2 atoms)
	// 	arma::field<arma::cx_cube> psi; // electronic wave functions corresponding to electronic states using the reduced graphene unit cell (2 atoms)
	// }

	// struct to bundle data and metadata of coulomg interaction fourier transform (vq)
	struct vq_struct
	{
		arma::cx_cube data; // actual data of vq in the format of (iq,mu,atom_pair_index) where atom pair index is aa=0, ab=1, ba=2, bb=3
		std::array<int,2> iq_range; // range of iq values in the half-open range format [a,b)
		std::array<int,2> mu_range; // range of mu values in the half-open range format [a,b)
		int nq, n_mu; // number of iq and mu elements
	};
	// instantiation of vq_struct to hold data of vq calculated via calculate_vq function
	vq_struct _vq;

	// struct to bundle data and metadata of electronic state polarization (PI)
	struct PI_struct
	{
		arma::mat data; // actual data of PI in the format of (iq,mu)
		std::array<int,2> iq_range; // range of iq values in the half-open range format [a,b)
		std::array<int,2> mu_range; // range of mu values in the half-open range format [a,b)
		int nq, n_mu; // number of iq and mu elements
	};
	// instantiation of PI_struct to hold data of PI calculated via calculate_polarization function
	PI_struct _PI;

	// struct to bundle data and metadata of dielectric function (epsilon)
	struct epsilon_struct
	{
		arma::mat data; // actual data of dielectric function in the format of (iq,mu)
		std::array<int,2> iq_range; // range of iq values in the half-open range format [a,b)
		std::array<int,2> mu_range; // range of mu values in the half-open range format [a,b)
		int nq, n_mu; // number of iq and mu elements
	};
	// instantiation of epsilon_struct to hold data of dielectric function calculated via calculate_dielectric function
	epsilon_struct _eps;

	int _i_sub = 0; // index of the selected subband from _valleys_K2 vector

	arma::cx_vec _epsilon; // static dielectric function


public:
	//constructor
	cnt(){};
	
	// set the output directory and the output file name
	void process_command_line_args(int argc, char* argv[]);

	// calculate the parameters of the cnt
	void get_parameters();
	
	// calculates position of atoms and reciprocal lattice vectors
	void get_atom_coordinates();
	
	// calculate electron dispersion energies using full unit cell (2*Nu atoms)
	void electron_full();
	
	// calculate electron dispersion energies using the reduced graphene unit cell (2 atoms)
	void electron_K1_extended();

	// calculate electron dispersion energies using the K2-extended representation
	void electron_K2_extended();

	// find valley ik and i_mu indices in K2-extended representation
	void find_K2_extended_valleys();

	// find ik values that are energetically relevant around the bottom of the valley
	void find_relev_ik_range(double delta_energy);

	// get ik of valence band state by taking care of wrapping around K2-extended zone
	int get_ikv(int ik_c, int ik_cm)
	{
		int ik_v = ik_c - ik_cm;
		while (ik_v >= _ik_max_K2){
			ik_v -= _nk_K2;
		}
		while (ik_v < _ik_min_K2){
			ik_v += _nk_K2;
		}
		return ik_v;
	};

	// get ik of conduction band state by taking care of wrapping around K2-extended zone
	int get_ikc(int ik_v, int ik_cm)
	{
		int ik_c = ik_v + ik_cm;
		while (ik_c >= _ik_max_K2){
			ik_c -= _nk_K2;
		}
		while (ik_c < _ik_min_K2){
			ik_c += _nk_K2;
		}
		return ik_c;
	};

	// helper function to check if a number is inside another range
	bool in_range(const int& guest, const std::array<int,2>& host) const
	{
		if (guest < host[0]){
			return false;
		}
		if (guest >= host[1]){
			return false;
		}
		return true;
	};
	// helper function to check if a range is inside another range
	bool in_range(const std::array<int,2>& guest, const std::array<int,2>& host) const
	{
		if (guest[0] < host[0]){
			return false;
		}
		if (guest[0] >= host[1]){
			return false;
		}
		if (guest[1] < host[0]){
			return false;
		}
		if (guest[1] > host[1]){
			return false;
		}
		return true;
	};
	
	// fourier transformation of the coulomb interaction a.k.a v(q)
	vq_struct calculate_vq(const std::array<int,2> iq_range, const std::array<int,2> mu_range, const unsigned int no_of_cnt_unit_cells);

	// polarization of electronic states a.k.a PI(q)
	PI_struct calculate_polarization(const std::array<int,2> iq_range, const std::array<int,2> mu_range);

	// dielectric function a.k.a eps(q)
	epsilon_struct calculate_dielectric(const std::array<int,2> iq_range, const std::array<int,2> mu_range);

	// call this to do all the calculations at once
	void calculate_exciton_dispersion();
};

#endif // end _cnt_h_
