#ifndef __CNT_H__
#define __CNT_H__

#include <stdio.h>
#include <string>
#include <math.h>
#include <complex>

#include "constants.h"
#include "nr3.h"

class cnt
{
private:

	// const double a_cc = pow(1.42,-10); // carbon-carbon distance [meters]
	const double a_cc = 1.42e-10; // carbon-carbon distance [meters]
	const double a_l = sqrt(3.0)*a_cc; // graphene lattice constants [meters]
	
	const double e2p = 0.0; // tight binding constants
	const double t0 = 2.7 * constants::eV; // tight binding constants
	const double s0 = 0.0; // tight binding constants

	const double Upp = 11.3 * constants::eV; // constant used in the Ohno potential

	int n, m; // chirailty parameters
	int length_in_cnt_unit_cell; // length of cnt in units of cnt unit cell.
	std::string name; //cnt name

	nr::vec_doub a1, a2; // real space lattice primitive vectors
	nr::vec_doub b1, b2; // reciprocal lattice primitive vectors
	nr::vec_doub aCC_vec; // vector between two neighboring carbon atoms
	nr::vec_doub ch_vec; // chirality vector

	double ch_len; // length of the chirality vector
	double radius; // radius of the cnt

	nr::vec_doub t_vec; // translational vector for cnt unit cell in unrolled graphene sheet
	nr::vec_doub t_vec_3d; // translational vector for cnt unit cell in rolled graphene sheed (3d)

	int Nu; // number of graphene unit cells in cnt unit cell.

	nr::vec_doub K1; // cnt reciprocal lattice vector in the circumferencial direction
	nr::vec_doub K2; // cnt reciprocal lattice vector along the cnt axis
	nr::vec_doub dk_l; // delta_k in the longitudinal direction with respect to cnt axis

	nr::mat_doub pos_a, pos_b; // position of atoms in A and B sites
	nr::mat_doub pos_2d, pos_3d; // position of all atoms in cnt unit cell in 2d and in 3d space
	// nr::mat_doub pos_aa, pos_ab, pos_ba, pos_bb; // distance between atoms and A and B sites which conserves lattice symmetry.

	nr::mat_doub ek;
	// nr::mat_complex psi;


public:
	cnt(const std::string &in_name, const int in_n, const int in_m, const int in_length);

	void geometry(); // calculates position of atoms
	void electron(); // electron dispersion energies
};

#endif
