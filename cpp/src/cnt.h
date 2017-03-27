#ifndef __CNT_H__
#define __CNT_H__

#include <stdio.h>
#include <string>
#include <math.h>
#include <complex>
#include <valarray>

#include <boost/numeric/ublas/vector.hpp>	// for boost::numeric::ublas vector definition
#include <boost/numeric/ublas/matrix.hpp>	// for boost::numeric::ublas matrix definition

#include "constants.h"

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

	boost::numeric::ublas::vector<double> a1, a2; // real space lattice primitive vectors
	boost::numeric::ublas::vector<double> b1, b2; // reciprocal lattice primitive vectors
	boost::numeric::ublas::vector<double> aCC_vec; // vector between two neighboring carbon atoms
	boost::numeric::ublas::vector<double> ch_vec; // chirality vector

	double ch_len; // length of the chirality vector
	double radius; // radius of the cnt

	boost::numeric::ublas::vector<double> t_vec; // translational vector for cnt unit cell in unrolled graphene sheet
	boost::numeric::ublas::vector<double> t_vec_3d; // translational vector for cnt unit cell in rolled graphene sheed (3d)

	int Nu; // number of graphene unit cells in cnt unit cell.

	boost::numeric::ublas::vector<double> K1; // cnt reciprocal lattice vector in the circumferencial direction
	boost::numeric::ublas::vector<double> K2; // cnt reciprocal lattice vector along the cnt axis
	boost::numeric::ublas::vector<double> dk_l; // delta_k in the longitudinal direction with respect to cnt axis

	boost::numeric::ublas::matrix<double> pos_a, pos_b; // position of atoms in A and B sites
	boost::numeric::ublas::matrix<double> pos_2d, pos_3d; // position of all atoms in cnt unit cell in 2d and in 3d space
	// boost::numeric::ublas::matrix<double> pos_aa, pos_ab, pos_ba, pos_bb; // distance between atoms and A and B sites which conserves lattice symmetry.

	boost::numeric::ublas::matrix<double> ek;
	boost::numeric::ublas::matrix<std::complex<double>> psi;


public:
	cnt();
	cnt(const char *in_name, int in_n, int in_m, int in_length);

	void geometry(); // calculates position of atoms
	void electron(); // electron dispersion energies
};

#endif
