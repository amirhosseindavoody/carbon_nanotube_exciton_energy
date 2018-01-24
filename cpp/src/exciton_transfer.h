#ifndef _exciton_transfer_h_
#define _exciton_transfer_h_

#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>
#include <type_traits>

#include "cnt.h"

class exciton_transfer
{
private:
  typedef std::experimental::filesystem::directory_entry t_directory; // custom directory type
  t_directory _directory; // this is the address of the directory that the simulation data is stored

  double _temperature; // temperature of the system to calculate thermal distribution for excitons and other particles
  // double _c2c_distance; // center to center distance between the cnts.
  // double _theta; // angle between the axis of the cnts
  double _broadening_factor; // broadening factor used in the lorenzian to simulate dirac delta function
  // std::array<double,2> _shifts; // the amount that each cnt is shifted from the center.
  std::array<const cnt*,2> _cnts; // array of pointers to the target excitons

  // function to return the lorentzian based on the broadening factor
  const double lorentzian(const double& energy)
  {
    return constants::inv_pi*_broadening_factor/(energy*energy + _broadening_factor*_broadening_factor);
  };

  // prepare the output directory by deleting its previous content or creating it
  t_directory prepare_directory();

public:
  // constructor
  exciton_transfer(const cnt& cnt1, const cnt& cnt2)
  {
    _directory = prepare_directory();

    _cnts = {&cnt1, &cnt2};
    _temperature = 300;
    _broadening_factor = 4.e-3*constants::eV;

    std::cout << "\n...exciton transfer parameters:\n";
    std::cout << "temperature: " << _temperature << " [Kelvin]\n";
    std::cout << "energy broadening factor: " << _broadening_factor/constants::eV << " [eV]\n";

  };

  // calculate first order transfer rate
  double first_order(const double& z_shift, const std::array<double,2> axis_shifts, const double& theta);

  // calculate and plot Q matrix element between two exciton bands
  void save_Q_matrix_element(const int i_n_principal, const int f_n_principal);

  // calculate and plot J matrix element between two exciton bands
  void save_J_matrix_element(const int i_n_principal, const int f_n_principal);

  // struct to bundle information about the excitonic states that are relevant
  struct ex_state
  {
    ex_state(const cnt::exciton_struct& m_exciton, const int& m_ik_cm_idx, const int& m_i_principal):
      exciton(&m_exciton), cnt_obj(m_exciton.cnt_obj), elec_struct(m_exciton.elec_struct)
    {
      ik_cm = m_ik_cm_idx+m_exciton.ik_cm_range[0];
      i_principal = m_i_principal;
      ik_cm_idx = m_ik_cm_idx;
      energy = m_exciton.energy(ik_cm_idx,i_principal);
    };

    const cnt::exciton_struct* exciton; // reference to the exciton struct that owns the state
    const cnt* cnt_obj; // reference to the cnt object owning the exciton state
    const cnt::el_energy_struct* elec_struct; // reference to the el_energy_struct that is used to calculate the exciton dipersion
    int ik_cm_idx; // index of the ik_cm state in the exciton.energy and exciton.psi matrices
    int i_principal; // the principal quantum number of the state in the exciton.energy and exciton.psi matrices

    double energy; // energy of the exciton state
    int ik_cm; // value of the ik_cm for the state
    // const int mu_cm=0; // value of the mu for the state

    // access to the whole exciton state wavefunction
    arma::cx_vec psi() const
    {
      return exciton->psi.slice(ik_cm_idx).col(i_principal);
    };

    // access to individual elements of exciton state wavefunction
    std::complex<double> psi(int ik_c_idx) const
    {
      return exciton->psi(ik_c_idx,i_principal,ik_cm_idx);
    };

    const arma::umat& ik_idx() const
    {
      return exciton->ik_idx.slice(ik_cm_idx);
    }

    unsigned int ik_idx(const int& j, const int& i_n_principal) const
    {
      return exciton->ik_idx(j,i_n_principal, ik_cm_idx);
    }

    const arma::vec& dk_l() const
    {
      return *(exciton->dk_l);
    }

    const arma::vec K_cm() const
    {
      return ik_cm*(*(exciton->dk_l));
    } 

  };

  // get the energetically relevant states in the form a vector of ex_state structs
  std::vector<ex_state> get_relevant_states(const cnt::exciton_struct& exciton, const double min_energy);

  // struct to bundle information about initial and final states that match energetically
  struct matching_states
  {
    matching_states(const ex_state& d_state, const ex_state& a_state): i(d_state), f(a_state) 
    {};
    const ex_state i; // initial exciton state
    const ex_state f; // final exciton state
  };

  // calculate Q()
  std::complex<double> calculate_Q(const matching_states& pair) const;

  // calculate J()
  std::complex<double> calculate_J(const matching_states& pair, const std::array<double,2>& shifts_along_axis, const double& z_shift, const double& angle) const;

  // match states based on energies
  std::vector<matching_states> match_states(const std::vector<ex_state>& d_relevant_states, const std::vector<ex_state>& a_relevant_states)
  {
    std::vector<matching_states> matched;
    
    for (const auto& d_state: d_relevant_states)
    {
      for (const auto& a_state: a_relevant_states)
      {
        if (is_matched(d_state, a_state))
        {
          matched.emplace_back(matching_states(d_state,a_state));
        }
      }
    }



    std::cout << "\n...calculated pairs of donor and acceptor states\n";
    std::cout << "number of pairs: " << matched.size() << " out of " << a_relevant_states.size()*d_relevant_states.size() << "\n\n";
    // for (const auto& pair:matched)
    // { 
    //   double delta_e = pair.i.energy-pair.f.energy;
    //   std::cout << "delta energy:" << delta_e/constants::eV << " [eV]   , lorentzian:" << lorentzian(delta_e)/lorentzian(0) << "\n";
    // }
    
    return matched;
  };

  // function to check if a state is energetically matched using a lorentzian
  bool is_matched(const ex_state& d_state, const ex_state& a_state)
  {
    double delta_e = d_state.energy - a_state.energy;
    if (lorentzian(delta_e) > 1.e-2*lorentzian(0))
    {
      return true;
    }
    return false;
  };

  // match states based on energies
  std::vector<matching_states> match_all_states(const std::vector<ex_state>& d_relevant_states, const std::vector<ex_state>& a_relevant_states)
  {
    std::vector<matching_states> matched;
    
    for (const auto& d_state: d_relevant_states)
    {
      for (const auto& a_state: a_relevant_states)
      {
        matched.emplace_back(matching_states(d_state,a_state));
      }
    }



    std::cout << "\n...matched all possible pairs of donor and acceptor states\n";
    std::cout << "number of pairs: " << matched.size() << "\n\n";
    
    return matched;
  };

};

#endif //_exciton_transfer_h_