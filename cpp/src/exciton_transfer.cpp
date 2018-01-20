#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>
#include <stdexcept>

#include "exciton_transfer.h"
#include "cnt.h"
#include "constants.h"

// prepare the output directory by deleting its previous content or creating it
exciton_transfer::t_directory exciton_transfer::prepare_directory()
{
  namespace fs = std::experimental::filesystem;
  t_directory directory;

  directory.assign("/Users/amirhossein/research/exciton_transfer");
  std::cout << "output_directory: " << directory.path() << std::endl;

  if (not fs::exists(directory.path()))
  {
    std::cout << "warning: output directory does NOT exist!!!" << std::endl;
    std::cout << "output directory: " << directory.path() << std::endl;
    fs::create_directories(directory.path());
  }

  if (fs::is_directory(directory.path()))
  {
    if (not fs::is_empty(directory.path()))
    {
      std::cout << "warning: output directory is NOT empty!!!" << std::endl;
      std::cout << "output directory: " << directory.path() << std::endl;
      std::cout << "deleting the existing directory!!!" << std::endl;
      fs::remove_all(directory.path());
      fs::create_directories(directory.path());
    }
  }
  else
  {
    std::cout << "error: output path is NOT a directory!!!" << std::endl;
    std::cout << "output path: " << directory.path() << std::endl;
    throw std::invalid_argument("The input value for output directory is not acceptable.");
  }

  return directory;

};

// calculate and plot Q matrix element between two exciton bands
void exciton_transfer::save_Q_matrix_element(const int i_n_principal, const int f_n_principal)
{
  // set the donor and accdeptor cnts
  const cnt& init_cnt = *_cnts[0];
  const cnt& final_cnt = *_cnts[1];
  
  // set the donor and acceptor exciton structs
  const cnt::exciton_struct& i_exciton = init_cnt.A2_singlet();
  const cnt::exciton_struct& f_exciton = final_cnt.A2_singlet();

  std::cout << "\n...calculating Q matrix element\n";
  std::cout << "initial exciton ik_cm_range: [" << i_exciton.ik_cm_range[0] << "," << i_exciton.ik_cm_range[1] << "]\n";
  std::cout << "final exciton ik_cm_range: [" << f_exciton.ik_cm_range[0] << "," << f_exciton.ik_cm_range[1] << "]\n";

  // some error checking
  if ((i_n_principal >= i_exciton.n_principal) and (f_n_principal >= f_exciton.n_principal)){
    throw std::invalid_argument("i_n_principal or f_n_principla are two large.");
  }

  // find lists of relevant states in donor and acceptor excitons
  auto get_exciton_band = [](const auto& exciton, const int i_principal){
    std::vector<ex_state> states;
    for (int ik_cm_idx=0; ik_cm_idx<exciton.nk_cm; ik_cm_idx++)
    {
      ex_state state(exciton,ik_cm_idx,i_principal);
      states.emplace_back(state);
    }
    return states;
  };
  const auto i_relevant_states = get_exciton_band(i_exciton,i_n_principal);
  const auto f_relevant_states = get_exciton_band(f_exciton,f_n_principal);


  // match the states based on their energy
  std::vector<matching_states> state_pairs = match_states(i_relevant_states, f_relevant_states, true);

  // get min of ik_cm_idx from relevant states indices
  auto get_min_ik_cm_idx = [](const auto& relevant_states){
    int min_idx = 1.e9;
    for (const auto& state: relevant_states)
    {
      min_idx = (min_idx < state.ik_cm_idx) ? min_idx : state.ik_cm_idx;
    }
    return min_idx;
  };
  // get max of ik_cm_idx from relevant states indices
  auto get_max_ik_cm_idx = [](const auto& relevant_states){
    int max_idx = -1.e9;
    for (const auto& state: relevant_states)
    {
      max_idx = (max_idx > state.ik_cm_idx) ? max_idx : state.ik_cm_idx;
    }
    return max_idx;
  };
  int i_min_idx = get_min_ik_cm_idx(i_relevant_states);
  int i_max_idx = get_max_ik_cm_idx(i_relevant_states);
  int f_min_idx = get_min_ik_cm_idx(f_relevant_states);
  int f_max_idx = get_max_ik_cm_idx(f_relevant_states);

  arma::cx_mat Q_mat(i_max_idx-i_min_idx+1, f_max_idx-f_min_idx+1, arma::fill::zeros);
  arma::vec init_ik_cm(i_max_idx-i_min_idx+1, arma::fill::zeros);
  arma::vec final_ik_cm(f_max_idx-f_min_idx+1, arma::fill::zeros);
  for (const auto& pair:state_pairs)
  { 
    // std::cout << "pair.i.psi: " << arma::size(pair.i.exciton.psi.slice(pair.i.ik_cm_idx).col(pair.i.i_principal)) << "\n";
    // std::cout << "pair.i.psi: " << arma::size(pair.i.psi()) << "\n";
    Q_mat(pair.i.ik_cm_idx-i_min_idx, pair.f.ik_cm_idx-f_min_idx) = calculate_Q(pair);
    init_ik_cm(pair.i.ik_cm_idx-i_min_idx) = pair.i.ik_cm;
    final_ik_cm(pair.f.ik_cm_idx-f_min_idx) = pair.f.ik_cm;
  }

  arma::mat tmp;
  

  // save the real part matrix element Q to a file
  tmp= arma::real(Q_mat);
  std::string filename = _directory.path() / "matrix_element_q.real.dat";
  tmp.save(filename,arma::arma_ascii);

  // save the imaginary part matrix element Q to a file
  tmp= arma::imag(Q_mat);
  filename = _directory.path() / "matrix_element_q.imag.dat";
  tmp.save(filename,arma::arma_ascii);

  // save ik_cm of the inital states
  filename = _directory.path() / "matrix_element_q.init_ik_cm.dat";
  init_ik_cm.save(filename,arma::arma_ascii);

  // save ik_cm of the final states
  filename = _directory.path() / "matrix_element_q.final_ik_cm.dat";
  final_ik_cm.save(filename,arma::arma_ascii);

  std::cout << "\n...calculated and saved Q matrix element\n";

};

// calculate Q()
std::complex<double> exciton_transfer::calculate_Q(const matching_states& pair)
{
  // lambda function to calculate part of the Q matrix element that relates to each pair.
  auto Q_partial = [](const ex_state& state)
  {
    const int ic = 1;
    const int iv = 0;

    const int iA = 0;
    const int iB = 1;

    const arma::vec dA = {0,0};
    const arma::vec& dB = state.exciton.aCC_vec;
    const arma::cx_vec exp_factor({std::exp(std::complex<double>(0.,-1.)*arma::dot(state.ik_cm*state.exciton.dk_l,dA)),\
                                   std::exp(std::complex<double>(0.,-1.)*arma::dot(state.ik_cm*state.exciton.dk_l,dB))});

    std::complex<double> Q_partial = 0;
    for (int ik_c_idx=0; ik_c_idx<state.exciton.nk_c; ik_c_idx++)
    {
      const arma::cx_vec& Cc = state.elec_struct.wavefunc(state.ik_idx(1,ik_c_idx)).slice(state.ik_idx(0,ik_c_idx)).col(ic);
      const arma::cx_vec& Cv = state.elec_struct.wavefunc(state.ik_idx(3,ik_c_idx)).slice(state.ik_idx(2,ik_c_idx)).col(iv);
      Q_partial += state.psi(ik_c_idx)*arma::accu(Cc%arma::conj(Cv)%exp_factor);
    }

    return Q_partial;

  };

  return std::conj(Q_partial(pair.i))*Q_partial(pair.f);

};