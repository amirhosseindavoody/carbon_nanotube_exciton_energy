#ifndef _exciton_transfer_h_
#define _exciton_transfer_h_

#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>

#include "cnt.h"

class exciton_transfer
{
private:  
  std::experimental::filesystem::directory_entry _directory; // this is the address of the directory that the simulation data is stored

  double _temperature; // temperature of the system to calculate thermal distribution for excitons and other particles
  double _c2c_distance; // center to center distance between the cnts.
  double _theta; // angle between the axis of the cnts
  double _broadening_factor; // broadening factor used in the lorenzian to simulate dirac delta function
  std::array<double,2> _shifts; // the amount that each cnt is shifted from the center.
  std::array<const cnt*,2> _cnts; // array of pointers to the target excitons

  // function to return the lorentzian based on the broadening factor
  const double lorentzian(const double& energy)
  {
    return constants::inv_pi*_broadening_factor/(energy*energy + _broadening_factor*_broadening_factor);
  }

  // prepare the output directory _directory
  void prepare_directory();

public:
  // constructor
  exciton_transfer(const cnt& cnt1, const cnt& cnt2)
  {
    prepare_directory();

    _cnts = {&cnt1, &cnt2};
    _temperature = 300;
    _broadening_factor = 4.e-3*constants::eV;

    std::cout << "\n...exciton transfer parameters:\n";
    std::cout << "temperature: " << _temperature << " [Kelvin]\n";
    std::cout << "energy broadening factor: " << _broadening_factor/constants::eV << " [eV]\n";

  };

 // calculate first order transfer rate
void first_order()
{
  _c2c_distance = 1.5e-9;
  _theta = constants::pi/2.;
  _shifts = {0., 0.};

  std::cout << "\ncenter to center distance: " << _c2c_distance*1.e9 << " [nm]\n";
  std::cout << "cnt1 radius: " << _cnts[0]->radius()*1.e9 << " [nm], cnt2 radius: " << _cnts[1]->radius()*1.e9 << " [nm]\n";
  std::cout << "wall to wall distance: " << (_c2c_distance - _cnts[0]->radius() - _cnts[1]->radius())*1.e9 << " [nm]\n";

  // calculate the crossing points and same-energy points with the same energy between cnt_1 and cnt_2

  const cnt& donor = *_cnts[0];
  const cnt& acceptor = *_cnts[1];

  double min_energy = donor.A2_singlet().energy.min();
  std::cout << "\nminimum exciton energy in donor cnt is: " << min_energy/constants::eV << " [eV]\n";

  const cnt::exciton_struct& d_exciton = donor.A2_singlet();
  const cnt::exciton_struct& a_exciton = acceptor.A2_singlet();

  // find lists of relevant states in donor and acceptor excitons
  std::vector<std::array<int,2>> d_relevant_states_indices = get_relevant_states(d_exciton,min_energy);
  std::vector<std::array<int,2>> a_relevant_states_indices = get_relevant_states(a_exciton,min_energy);

  // match the states based on their energy
  std::vector<matching_states> state_pairs = match_states(d_exciton, a_exciton, d_relevant_states_indices, a_relevant_states_indices);

  auto get_min_idx = [](const auto& relevant_state_indices){
    int min_idx = 1.e9;
    for (const auto& idx: relevant_state_indices)
    {
      min_idx = (min_idx < idx[0]) ? min_idx : idx[0];
    }
    return min_idx;
  };
  
  auto get_max_idx = [](const auto& relevant_state_indices){
    int max_idx = -1.e9;
    for (const auto& idx: relevant_state_indices)
    {
      max_idx = (max_idx > idx[0]) ? max_idx : idx[0];
    }
    return max_idx;
  };

  int d_min_idx = get_min_idx(d_relevant_states_indices);
  int d_max_idx = get_max_idx(d_relevant_states_indices);
  int a_min_idx = get_min_idx(a_relevant_states_indices);
  int a_max_idx = get_max_idx(a_relevant_states_indices);

  arma::cx_mat Q_mat(d_max_idx-d_min_idx+1, a_max_idx-a_min_idx+1, arma::fill::zeros);

  for (const auto& pair:state_pairs)
  { 
    std::complex<double> Q = calculate_Q(pair);
    // Q_mat(pair.i_state_idx[0]-d_min_idx, pair.f_state_idx[0]-a_min_idx) = Q;
    Q_mat(pair.i_state_idx[0]-d_min_idx, pair.f_state_idx[0]-a_min_idx) = Q;
  }

  arma::mat tmp_Q(arma::abs(Q_mat));

  // save the matrix element Q to a file
  std::string filename = _directory.path() / "matrix_element_q.dat";
  tmp_Q.save(filename,arma::arma_ascii);

  std::cout << "\n...calculated Q matrix element\n";

};

// calculate and plot Q matrix element between two exciton bands
void save_Q_matrix_element(const int i_n_principal, const int f_n_principal);

// get the index of relevant energy states in the form [[ik_cm_idx, i_n], ... ]
std::vector<std::array<int,2>> get_relevant_states(const cnt::exciton_struct& exciton, const double min_energy)
{
  const double threshold_population = 1.e-3;
  const double threshold_energy = min_energy+std::abs(std::log(threshold_population) * constants::kb*_temperature);

  std::vector<std::array<int,2>> relevant_state_indices;
  
  for (int i_n=0; i_n<exciton.n_principal; i_n++)
  {
    for (int ik_cm_idx=0; ik_cm_idx<exciton.nk_cm; ik_cm_idx++)
    {
      if (exciton.energy(ik_cm_idx,i_n) <= threshold_energy)
      {
        relevant_state_indices.push_back({ik_cm_idx, i_n});
      }
    }
  }

  // sort relevant states in order of their energies
  std::sort(relevant_state_indices.begin(),relevant_state_indices.end(), \
              [&](const auto& s1, const auto& s2) {
                return exciton.energy(s1[0],s1[1]) < exciton.energy(s2[0],s2[1]);
              }
            );

  // calculate normalization factor
  double normalization_factor = 0.;
  for (const auto& idx: relevant_state_indices)
  {
    double delta_e = (exciton.energy(idx[0],idx[1])-min_energy);
    normalization_factor += std::exp(-delta_e/(constants::kb*_temperature));
  }


  std::cout << "\n...calculated relevant states\n";
  std::cout << "number of relevant states: " << relevant_state_indices.size() << "\n";
  // for (const auto& idx: relevant_state_indices)
  // {
  //   double delta_e = (exciton.energy(idx[0],idx[1])-min_energy);
  //   std::cout << "[" << idx[0] << "," << idx[1] <<"] --> energy:" << delta_e/constants::eV \
  //             << " , population:" << std::exp(-delta_e/(constants::kb*_temperature))/normalization_factor << "\n";
  // }
  
  return relevant_state_indices;
};

// struct to bundle information about initial and final states that match energetically
struct matching_states
{
  matching_states(const cnt::exciton_struct& d_exciton, const cnt::exciton_struct& a_exciton):
    i_exciton(d_exciton), f_exciton(a_exciton) {};
  std::array<int,2> i_state_idx, f_state_idx; // ik_cm_idx and mu_cm_idx of the initial and final states in the exciton_struct
  double i_energy, f_energy; // energy of the initial and final states
  const cnt::exciton_struct &i_exciton, &f_exciton; // references to the initial and final excitons so that we don't have to manually carry it everywhere
  int i_ik_cm=0, f_ik_cm=0; // ik_cm of the initial and final states

  // function to set properties of initial state of the pair
  void set_i_state(const std::array<int,2>& d_state)
  {
    i_state_idx = d_state;
    i_energy = i_exciton.energy(i_state_idx[0],i_state_idx[1]);
    i_ik_cm = i_state_idx[0]+i_exciton.ik_cm_range[0];
    // std::cout << "i_ik_cm:" << i_ik_cm << ", i_state_idx[0]:" << i_state_idx[0] << ", i_exciton.ik_cm_range[0]:" << i_exciton.ik_cm_range[0] << "\n";
  }

  // function to set properties of final state of the pair
  void set_f_state(const std::array<int,2>& a_state)
  {
    f_state_idx = a_state;
    f_energy = f_exciton.energy(f_state_idx[0],f_state_idx[1]);
    f_ik_cm = f_state_idx[0]+f_exciton.ik_cm_range[0];
  }
};

// calculate Q()
std::complex<double> calculate_Q(const matching_states& pair)
{
  const int ic = 1;
  const int iv = 0;

  const int iA = 0;
  const int iB = 1;


  // calculate Q for initial state
  const std::array<int,2>& i_state_idx = pair.i_state_idx;
  const cnt::exciton_struct& i_exciton = pair.i_exciton;
  const arma::cx_vec& Ai = i_exciton.psi.slice(i_state_idx[0]).col(i_state_idx[1]);
  const arma::umat& i_ik_idx = i_exciton.ik_idx.slice(i_state_idx[0]); // (j, i_elec_state)

  const arma::vec dA = {0,0};
  const arma::vec& dB = i_exciton.aCC_vec;
  const arma::cx_vec i_exp_factor({std::exp(std::complex<double>(0.,-1.)*arma::dot(pair.i_ik_cm*i_exciton.dk_l,dA)),\
                                   std::exp(std::complex<double>(0.,-1.)*arma::dot(pair.i_ik_cm*i_exciton.dk_l,dB))});
  
  std::complex<double> Q_i = 0;
  for (int ik_c_idx=0; ik_c_idx<i_exciton.nk_c; ik_c_idx++)
  {
    const arma::cx_vec& Cc = i_exciton.elec_struct.wavefunc(i_ik_idx(1,ik_c_idx)).slice(i_ik_idx(0,ik_c_idx)).col(ic);
    const arma::cx_vec& Cv = i_exciton.elec_struct.wavefunc(i_ik_idx(3,ik_c_idx)).slice(i_ik_idx(2,ik_c_idx)).col(iv);
    Q_i += Ai(ik_c_idx)*arma::accu(Cc%arma::conj(Cv)%i_exp_factor);
  }

  // calculate Q for the final state
  const std::array<int,2>& f_state_idx = pair.f_state_idx;
  const cnt::exciton_struct& f_exciton = pair.f_exciton;
  const arma::cx_vec& Af = f_exciton.psi.slice(f_state_idx[0]).col(f_state_idx[1]);
  const arma::umat& f_ik_idx = f_exciton.ik_idx.slice(f_state_idx[0]); // (j, i_elec_state)
  
  const arma::vec dAp = {0,0};
  const arma::vec& dBp = f_exciton.aCC_vec;
  const arma::cx_vec f_exp_factor({std::exp(std::complex<double>(0.,-1.)*arma::dot(pair.f_ik_cm*f_exciton.dk_l,dAp)),\
                                   std::exp(std::complex<double>(0.,-1.)*arma::dot(pair.f_ik_cm*f_exciton.dk_l,dBp))});
  std::complex<double> Q_f = 0;  
  for (int ik_cp_idx=0; ik_cp_idx<f_exciton.nk_c; ik_cp_idx++)
  {
    const arma::cx_vec& Ccp = f_exciton.elec_struct.wavefunc(f_ik_idx(1,ik_cp_idx)).slice(f_ik_idx(0,ik_cp_idx)).col(ic);
    const arma::cx_vec& Cvp = f_exciton.elec_struct.wavefunc(f_ik_idx(3,ik_cp_idx)).slice(f_ik_idx(2,ik_cp_idx)).col(iv);
    Q_f += Af(ik_cp_idx)*arma::accu(Ccp%arma::conj(Cvp)%arma::conj(f_exp_factor));
  }

  return std::conj(Q_i)*Q_f;
};

// match states based on energies
std::vector<matching_states> match_states(const cnt::exciton_struct& d_exciton, \
                                          const cnt::exciton_struct& a_exciton, \
                                          const std::vector<std::array<int,2>>& d_relevant_states_indices, \
                                          const std::vector<std::array<int,2>>& a_relevant_states_indices,
                                          const bool& get_all_possible=false)
{
  std::vector<matching_states> matched;
  
  for (const auto& d_state: d_relevant_states_indices)
  {
    for (const auto& a_state: a_relevant_states_indices)
    {
      if (is_matched(d_exciton, a_exciton, d_state, a_state) or get_all_possible)
      {
        matching_states my_pair(d_exciton, a_exciton);
        my_pair.set_i_state(d_state);
        my_pair.set_f_state(a_state);
        matched.push_back(my_pair);
      }
    }
  }



  std::cout << "\n...calculated pairs of donor and acceptor states\n";
  std::cout << "number of pairs: " << matched.size() << " out of " << a_relevant_states_indices.size()*d_relevant_states_indices.size() << "\n\n";
  // for (const auto& pair:matched)
  // { 
  //   double delta_e = pair.i_energy-pair.f_energy;
  //   std::cout << "delta energy:" << delta_e/constants::eV << " [eV]   , lorentzian:" << lorentzian(delta_e)/lorentzian(0) << "\n";
  // }
  
  return matched;
};

// function to check if a state is energetically matched using a lorentzian
bool is_matched(const cnt::exciton_struct& d_exciton, const cnt::exciton_struct& a_exciton, \
             const std::array<int,2>& d_state, const std::array<int,2>& a_state)
{
  double delta_e = d_exciton.energy(d_state[0],d_state[1]) - a_exciton.energy(a_state[0],a_state[1]);
  if (lorentzian(delta_e) > 1.e-2*lorentzian(0))
  {
    return true;
  }
  return false;
};

};

#endif //_exciton_transfer_h_