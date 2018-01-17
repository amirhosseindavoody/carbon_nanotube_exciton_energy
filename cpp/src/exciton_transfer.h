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
  double _temperature; // temperature of the system to calculate thermal distribution for excitons and other particles
  double _c2c_distance; // center to center distance between the cnts.
  double _theta; // angle between the axis of the cnts
  double _broadening_factor; // broadening factor used in the lorenzian to simulate dirac delta function
  std::array<double,2> _shifts; // the amount that each cnt is shifted from the center.
  std::array<const cnt*,2> _cnts; // array of pointers to the target excitons

  // find the points in dispersions of cnts that cross each other and conserve both energy and momentum
  void find_crossings();

  // function to return the lorentzian based on the broadening factor
  const double lorentzian(const double& energy)
  {
    return constants::inv_pi*_broadening_factor/(energy*energy + _broadening_factor*_broadening_factor);
  }

public:
  // constructor
  exciton_transfer(const cnt& cnt1, const cnt& cnt2)
  {
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

  for (const auto& pair:state_pairs)
  { 
    calculate_Q(pair);
  }

};

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
  // int ik_cm_i, ik_cm_f; // value of ik_cm for initial and final states
  // int mu_cm_i, mu_cm_f; // value of mu_cm for initial and final states

  // function to set properties of initial state of the pair
  void set_i_state(const std::array<int,2>& d_state)
  {
    i_state_idx = d_state;
    i_energy = i_exciton.energy(i_state_idx[0],i_state_idx[1]);
  }

  // function to set properties of final state of the pair
  void set_f_state(const std::array<int,2>& a_state)
  {
    f_state_idx = a_state;
    f_energy = f_exciton.energy(f_state_idx[0],f_state_idx[1]);
  }
};

// calculate Q()
void calculate_Q(const matching_states& pair)
{
  std::complex<double> Q = 0;
  // for (int ik_c=pair.; ik_c<)

  // std::cout << "\n...calculating Q matrix element\n";
}

// match states based on energies
std::vector<matching_states> match_states(const cnt::exciton_struct& d_exciton, \
                                          const cnt::exciton_struct& a_exciton, \
                                          const std::vector<std::array<int,2>>& d_relevant_states_indices, \
                                          const std::vector<std::array<int,2>>& a_relevant_states_indices )
{
  std::vector<matching_states> matched;
  
  for (const auto& d_state: d_relevant_states_indices)
  {
    for (const auto& a_state: a_relevant_states_indices)
    {
      if (is_matched(d_exciton, a_exciton, d_state, a_state))
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