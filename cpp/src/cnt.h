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

  std::vector<std::array<std::array<unsigned int, 2>, 2>> _valleys_K2; // index of valleys in K2-extended representation
  std::vector<std::vector<std::array<int,2>>> _relev_ik_range; // ik of relevant states in the following form [[[ik,mu],...], [[ik,mu],...]]

public:
  // struct to bundle information of electronic energy state
  struct el_energy_struct
  {
    std::string name; // human interpretable name describing the content of the struct
    int no_of_atoms; // number of atoms that are used in the wavefunction: it is 2 when graphen unit cell is used or 2*_Nu when full cnt unit cell is used.
    int no_of_bands; // number of bands for each choice of ik and mu: it is 2 when graphen unit cell is used and 2*_Nu when full cnt unit cell is used.
  	arma::cube energy; // energy of electronic states calculated using the reduced graphene unit cell (2 atoms)
  	arma::field<arma::cx_cube> wavefunc; // electronic wave functions corresponding to electronic states using the reduced graphene unit cell (2 atoms) \
                                            the wavefunc format is (wavefunc(mu-mu_range[0]))(iA,ic,ik-ik_range[0])
    std::array<int,2> ik_range;
    std::array<int,2> mu_range;
    int nk, n_mu; // number of elements in the range of ik and mu
  };
private:
  // instantiation of el_energy_struct within K2-extended representation
  el_energy_struct _elec_K2;

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
  // struct to bundle data and metadata of exciton
  struct exciton_struct
  {
    exciton_struct(const cnt* m_cnt): cnt_obj(*m_cnt), dk_l((m_cnt->_dk_l)), elec_struct(m_cnt->_elec_K2),
      aCC_vec(m_cnt->_aCC_vec) {};
    std::string name; // a human readable name for the exciton
    arma::mat energy; // exciton energy dispersion in the form (ik_cm, n) where n is the \
                         quantum number equivalent to principarl quantum number in hydrogen
    int spin; // exciton total spin
    std::array<int,2> ik_cm_range; // range of ik_cm values in the half-open range format [a,b)
    int mu_cm = 0; // exciton center of mass mu
    int n_principal=0; // number of states equivalent to the the principal quantum number in hydrogen atom
    int nk_c=0; // number of relevant states to make the exciton wave function
    int nk_cm=0; // number of ik_cm states
    arma::cx_cube psi; // exciton wavefunction in the form (ik_c_relev,n,ik_cm) therefore the first element is the \
                          weight of ik_c_relev state in the n-th eigen state with center-of-mass momentum ik_cm
    const el_energy_struct& elec_struct; // reference to the el_energy_struct that is used to calculate the exciton dipersion
    const arma::vec& dk_l; // reference to _dk_l vector in the owner cnt object
    const arma::vec& aCC_vec; // const reference to _aCC_vec in the owner cnt object

    arma::ucube ik_idx; // cube to hold index of kc and kv states for each element in psi.\
                           The cube has dimensions of (4, nk_relev, nk_cm) where \
                           element (j, i_elec_state, ik_cm_idx) shows index of ik_c, mu_c, ik_v, mu_v \
                           for the corresponding excitonic state: \
                           j=0 --> ik_c_idx, j=1 --> mu_c_idx, j=2 --> ik_v_idx, j=3 --> mu_v_idx
    
    const cnt& cnt_obj; // constant reference to the owner cnt object
    
    // assignment operator
    exciton_struct operator=(const exciton_struct& other)
    {
      exciton_struct tmp(&(other.cnt_obj));
      tmp.name = other.name;
      tmp.energy = other.energy;
      tmp.spin = other.spin;
      tmp.mu_cm = other.mu_cm;
      tmp.n_principal = other.n_principal;
      tmp.nk_c = other.nk_c;
      tmp.nk_cm = other.nk_cm;
      tmp.psi = other.psi;
      tmp.ik_idx = other.ik_idx;
      tmp.ik_cm_range = other.ik_cm_range;
      return tmp;
    };

    // copy constructor
    exciton_struct (const exciton_struct& other):cnt_obj(other.cnt_obj), dk_l((other.dk_l)), elec_struct(other.elec_struct),
      aCC_vec(other.aCC_vec)
    {
      (*this).name = other.name;
      (*this).energy = other.energy;
      (*this).spin = other.spin;
      (*this).mu_cm = other.mu_cm;
      (*this).n_principal = other.n_principal;
      (*this).nk_c = other.nk_c;
      (*this).nk_cm = other.nk_cm;
      (*this).psi = other.psi;
      (*this).ik_idx = other.ik_idx;
      (*this).ik_cm_range = other.ik_cm_range;
    };

    // move constructor
    exciton_struct (const exciton_struct&& other):cnt_obj(other.cnt_obj), dk_l((other.dk_l)), elec_struct(other.elec_struct),
      aCC_vec(other.aCC_vec)
    {
      (*this).name = other.name;
      (*this).energy = other.energy;
      (*this).spin = other.spin;
      (*this).mu_cm = other.mu_cm;
      (*this).n_principal = other.n_principal;
      (*this).nk_c = other.nk_c;
      (*this).nk_cm = other.nk_cm;
      (*this).psi = other.psi;
      (*this).ik_idx = other.ik_idx;
      (*this).ik_cm_range = other.ik_cm_range;
    };
  };

private:
  // vector of exciton structs to hold data of A and E type excitons
  std::vector<exciton_struct>  _excitons;

public:
  //constructor
  cnt(){};
  
  // set the output directory and the output file name
  void process_command_line_args(int argc, char* argv[]);

  // calculate the parameters of the cnt
  void get_parameters();
  
  // calculates position of atoms and reciprocal lattice vectors
  void get_atom_coordinates();
  
  // calculate electron energy dispersions in the K1-extended representation using full unit cell (2*Nu atoms)
  void electron_full_unit_cell();

  // calculate electron dispersion energies using the K2-extended representation
  void electron_K2_extended();

  // calculate electron dispersion energies for an input range of ik and mu
  el_energy_struct electron_energy(const std::array<int,2>& ik_range, const std::array<int,2>& mu_range, const std::string& name);

  // find valley ik and i_mu indices in K2-extended representation
  void find_valleys(const el_energy_struct& elec_struct);

  // find ik values that are energetically relevant around the bottom of the valley
  void find_relev_ik_range(double delta_energy, const el_energy_struct& elec_struct);

  // fourier transformation of the coulomb interaction a.k.a v(q)
  vq_struct calculate_vq(const std::array<int,2> iq_range, const std::array<int,2> mu_range, const unsigned int no_of_cnt_unit_cells);

  // polarization of electronic states a.k.a PI(q)
  PI_struct calculate_polarization(const std::array<int,2> iq_range, const std::array<int,2> mu_range, const el_energy_struct& elec_struct);

  // dielectric function a.k.a eps(q)
  epsilon_struct calculate_dielectric(const std::array<int,2> iq_range, const std::array<int,2> mu_range);

  // calculate exciton dispersion
  std::vector<exciton_struct> calculate_A_excitons(const std::array<int,2> ik_cm_range, const el_energy_struct& elec_struct);

  // call this to do all the calculations at once
  void calculate_exciton_dispersion();

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

  // getter function to access cnt radius
  const double& radius() const
  {
    return _radius;
  };

  // getter function to access all excitons
  const std::vector<exciton_struct>& excitons() const
  {
    return _excitons;
  };

  // getter function to access A2 singlet exciton
  const cnt::exciton_struct& A2_singlet() const
  {
    for (const auto& exciton: _excitons)
    {
      if (exciton.name == "A2 singlet exciton")
      {
        return exciton;
      }
    }
    throw std::logic_error("could not find A2 singlet exciton.\nInvestigate.\n");
  };

  // getter function to access A2 triplet exciton
  const cnt::exciton_struct& A2_triplet() const
  {
    for (const auto& exciton: _excitons)
    {
      if (exciton.name == "A2 triplet exciton")
      {
        return exciton;
      }
    }
    throw std::logic_error("could not find A2 triplet exciton.\nInvestigate.\n");
  };

  // getter function to access A1 exciton
  const cnt::exciton_struct& A1() const
  {
    for (const auto& exciton: _excitons)
    {
      if (exciton.name == "A1 exciton")
      {
        return exciton;
      }
    }
    throw std::logic_error("could not find A1 exciton.\nInvestigate.\n");
  };

  // helper class to monitor progress of the loops
  class progress_bar
  {
    private:
    float counter = 0;

    public:
    void step(const int& i, const int& i_max, const std::string& message, const int& increments=1)
    {
      if ((float(i)/float(i_max))>(counter/100.))
      {
        std::cout << message << ": " << int(counter) << "%\n";
        counter += increments;
      }
    };

    void reset()
    {
      counter = 0;
    }
  };
};

#endif // end _cnt_h_
