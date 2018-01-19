#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>
#include <stdexcept>

#include "exciton_transfer.h"
#include "cnt.h"
#include "constants.h"

// set the output directory and the output file name
void exciton_transfer::prepare_directory()
{
  namespace fs = std::experimental::filesystem;

  _directory.assign("/Users/amirhossein/research/exciton_transfer");
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
    throw std::invalid_argument("The input value for output directory is not acceptable.");
  }

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
    std::vector<std::array<int,2>> states_indices;
    for (int ik_cm_idx=0; ik_cm_idx<exciton.nk_cm; ik_cm_idx++)
    {
      states_indices.push_back({ik_cm_idx, i_principal});
    }
    return states_indices;
  };
  const auto i_relevant_states_indices = get_exciton_band(i_exciton,i_n_principal);
  const auto f_relevant_states_indices = get_exciton_band(f_exciton,f_n_principal);


  // match the states based on their energy
  std::vector<matching_states> state_pairs = match_states(i_exciton, f_exciton, i_relevant_states_indices, f_relevant_states_indices, true);

  // get min of ik_cm_idx from relevant states indices
  auto get_min_idx = [](const auto& relevant_state_indices){
    int min_idx = 1.e9;
    for (const auto& idx: relevant_state_indices)
    {
      min_idx = (min_idx < idx[0]) ? min_idx : idx[0];
    }
    return min_idx;
  };
  // get max of ik_cm_idx from relevant states indices
  auto get_max_idx = [](const auto& relevant_state_indices){
    int max_idx = -1.e9;
    for (const auto& idx: relevant_state_indices)
    {
      max_idx = (max_idx > idx[0]) ? max_idx : idx[0];
    }
    return max_idx;
  };
  int i_min_idx = get_min_idx(i_relevant_states_indices);
  int i_max_idx = get_max_idx(i_relevant_states_indices);
  int f_min_idx = get_min_idx(f_relevant_states_indices);
  int f_max_idx = get_max_idx(f_relevant_states_indices);

  arma::cx_mat Q_mat(i_max_idx-i_min_idx+1, f_max_idx-f_min_idx+1, arma::fill::zeros);
  arma::vec init_ik_cm(i_max_idx-i_min_idx+1, arma::fill::zeros);
  arma::vec final_ik_cm(f_max_idx-f_min_idx+1, arma::fill::zeros);
  for (const auto& pair:state_pairs)
  { 
    Q_mat(pair.i_state_idx[0]-i_min_idx, pair.f_state_idx[0]-f_min_idx) = calculate_Q(pair);
    init_ik_cm(pair.i_state_idx[0]-i_min_idx) = pair.i_ik_cm;
    final_ik_cm(pair.f_state_idx[0]-f_min_idx) = pair.f_ik_cm;
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