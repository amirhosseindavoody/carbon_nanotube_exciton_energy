#include <iostream>
#include <string>
#include <experimental/filesystem>
#include <armadillo>
#include <stdexcept>

#include "exciton_transfer.h"
#include "cnt.h"
#include "constants.h"
#include "progress.hpp"

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

}

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
  std::vector<matching_states> state_pairs = match_all_states(i_relevant_states, f_relevant_states);

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

  progress_bar prog(state_pairs.size(), "calculate Q");

  for (const auto& pair:state_pairs)
  { 
    Q_mat(pair.i.ik_cm_idx-i_min_idx, pair.f.ik_cm_idx-f_min_idx) = calculate_Q(pair);
    init_ik_cm(pair.i.ik_cm_idx-i_min_idx) = pair.i.ik_cm;
    final_ik_cm(pair.f.ik_cm_idx-f_min_idx) = pair.f.ik_cm;
    prog.step();
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

}

// calculate and plot J matrix element between two exciton bands
void exciton_transfer::save_J_matrix_element(const int i_n_principal, const int f_n_principal)
{
  // set geometrical properties
  double z_shift = 1.5e-9;
  double theta = 0;
  std::array<double,2> axis_shifts = {0., 0.};

  // set the donor and accdeptor cnts
  const cnt& init_cnt = *_cnts[0];
  const cnt& final_cnt = *_cnts[1];
  
  // set the donor and acceptor exciton structs
  const cnt::exciton_struct& i_exciton = init_cnt.A2_singlet();
  const cnt::exciton_struct& f_exciton = final_cnt.A2_singlet();

  std::cout << "\n...calculating J matrix element\n";
  // std::cout << "initial exciton ik_cm_range: [" << i_exciton.ik_cm_range[0] << "," << i_exciton.ik_cm_range[1] << "]\n";
  // std::cout << "final exciton ik_cm_range: [" << f_exciton.ik_cm_range[0] << "," << f_exciton.ik_cm_range[1] << "]\n";

  // some error checking
  if ((i_n_principal >= i_exciton.n_principal) and (f_n_principal >= f_exciton.n_principal)){
    throw std::invalid_argument("i_n_principal or f_n_principla are two large.");
  }

  // find lists of relevant states in donor and acceptor excitons
  auto get_exciton_band = [](const auto& exciton, const int i_principal){
    std::vector<ex_state> states;
    // int ik_cm_min = exciton.ik_cm_range[0];
    // int ik_cm_max = exciton.ik_cm_range[1];
    int ik_cm_min = -50;
    int ik_cm_max = -ik_cm_min+1;
    for (int ik_cm=ik_cm_min; ik_cm<ik_cm_max; ik_cm++)
    {
      ex_state state(exciton,ik_cm-exciton.ik_cm_range[0],i_principal);
      states.emplace_back(state);
    }
    return states;
  };
  const auto i_relevant_states = get_exciton_band(i_exciton,i_n_principal);
  const auto f_relevant_states = get_exciton_band(f_exciton,f_n_principal);

  std::cout << "initial exciton length: " << i_exciton.cnt_obj->length_in_meter()*1e9 << " [nm]\n";
  std::cout << "final exciton length: " << f_exciton.cnt_obj->length_in_meter()*1e9 << " [nm]\n";


  // match the states based on their energy
  std::vector<matching_states> state_pairs = match_all_states(i_relevant_states, f_relevant_states);

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

  arma::cx_mat J_mat(i_max_idx-i_min_idx+1, f_max_idx-f_min_idx+1, arma::fill::zeros);
  arma::vec init_ik_cm(i_max_idx-i_min_idx+1, arma::fill::zeros);
  arma::vec final_ik_cm(f_max_idx-f_min_idx+1, arma::fill::zeros);

  // int count = 0;
  // int number_of_pairs = state_pairs.size();
  progress_bar prog(state_pairs.size(),"calculate J");

  for (const auto& pair:state_pairs)
  { 
    // count++;
    // prog.step(count, number_of_pairs, "calculate J");
    prog.step();
    // std::cout << "calculate J: " << count << "/" << number_of_pairs << "\n";

    J_mat(pair.i.ik_cm_idx-i_min_idx, pair.f.ik_cm_idx-f_min_idx) = calculate_J(pair, axis_shifts, z_shift, theta);
    init_ik_cm(pair.i.ik_cm_idx-i_min_idx) = pair.i.ik_cm;
    final_ik_cm(pair.f.ik_cm_idx-f_min_idx) = pair.f.ik_cm;
  }

  arma::mat tmp;
  

  // save the real part matrix element Q to a file
  tmp= arma::real(J_mat);
  std::string filename = _directory.path() / "matrix_element_j.real.dat";
  tmp.save(filename,arma::arma_ascii);

  // save the imaginary part matrix element Q to a file
  tmp= arma::imag(J_mat);
  filename = _directory.path() / "matrix_element_j.imag.dat";
  tmp.save(filename,arma::arma_ascii);

  // save ik_cm of the inital states
  filename = _directory.path() / "matrix_element_j.init_ik_cm.dat";
  init_ik_cm.save(filename,arma::arma_ascii);

  // save ik_cm of the final states
  filename = _directory.path() / "matrix_element_j.final_ik_cm.dat";
  final_ik_cm.save(filename,arma::arma_ascii);

  std::cout << "\n...calculated and saved J matrix element\n";

}

// get the energetically relevant states in the form a vector of ex_state structs
std::vector<exciton_transfer::ex_state> exciton_transfer::get_relevant_states(const cnt::exciton_struct& exciton, const double min_energy)
{
  const double threshold_population = 1.e-3;
  const double threshold_energy = min_energy+std::abs(std::log(threshold_population) * constants::kb*_temperature);

  std::vector<ex_state> relevant_states;
  
  for (int i_n=0; i_n<exciton.n_principal; i_n++)
  {
    for (int ik_cm_idx=0; ik_cm_idx<exciton.nk_cm; ik_cm_idx++)
    {
      if (exciton.energy(ik_cm_idx,i_n) <= threshold_energy)
      {
        relevant_states.emplace_back(ex_state(exciton,ik_cm_idx,i_n));
      }
    }
  }

  // sort relevant states in order of their energies
  std::sort(relevant_states.begin(),relevant_states.end(), \
              [](const auto& s1, const auto& s2) {
                return s1.energy < s2.energy;
              }
            );

  // calculate normalization factor
  double normalization_factor = 0.;
  for (const auto& state: relevant_states)
  {
    double delta_e = (state.energy-min_energy);
    normalization_factor += std::exp(-delta_e/(constants::kb*_temperature));
  }


  // std::cout << "\n...calculated relevant states\n";
  // std::cout << "number of relevant states: " << relevant_states.size() << "\n";
  // for (const auto& state: relevant_states)
  // {
  //   double delta_e = (state.energy-min_energy);
  //   std::cout << "[" << state.ik_cm_idx << "," << state.i_principal <<"] --> energy:" << delta_e/constants::eV \
  //             << " , population:" << std::exp(-delta_e/(constants::kb*_temperature))/normalization_factor << "\n";
  // }
  
  return relevant_states;
}

// calculate Q()
std::complex<double> exciton_transfer::calculate_Q(const matching_states& pair) const
{
  // lambda function to calculate part of the Q matrix element that relates to each pair.
  auto Q_partial = [](const ex_state& state)
  {
    const int ic = 1;
    const int iv = 0;

    const arma::vec dA = {0,0};
    const arma::vec& dB = *(state.exciton->aCC_vec);
    const arma::cx_vec exp_factor({std::exp(std::complex<double>(0.,+1.)*arma::dot(state.ik_cm*state.dk_l(),dA)),\
                                   std::exp(std::complex<double>(0.,+1.)*arma::dot(state.ik_cm*state.dk_l(),dB))});

    std::complex<double> Q_partial = 0;
    for (int ik_c_idx=0; ik_c_idx<state.exciton->nk_c; ik_c_idx++)
    {
      const arma::cx_vec& Cc = state.elec_struct->wavefunc(state.ik_idx(1,ik_c_idx)).slice(state.ik_idx(0,ik_c_idx)).col(ic);
      const arma::cx_vec& Cv = state.elec_struct->wavefunc(state.ik_idx(3,ik_c_idx)).slice(state.ik_idx(2,ik_c_idx)).col(iv);
      Q_partial += state.psi(ik_c_idx)*arma::accu(Cc%arma::conj(Cv)%exp_factor);
    }

    return Q_partial;

  };

  double coeff = (std::pow(constants::q0,2)*pair.i.cnt_obj->Au()*pair.f.cnt_obj->Au())/\
                 (16*std::pow(constants::pi,3)*constants::eps0*pair.i.cnt_obj->radius()*pair.f.cnt_obj->radius()*\
                  std::sqrt(pair.i.cnt_obj->length_in_meter()*pair.f.cnt_obj->length_in_meter()));

  return std::complex<double>(coeff)*std::conj(Q_partial(pair.i))*Q_partial(pair.f);

}

// calculate J()
std::complex<double> exciton_transfer::calculate_J(const matching_states& pair, const std::array<double,2>& shifts_along_axis, const double& z_shift, const double& angle) const
{
  // make position of all atoms in the entire cnt length in 3d space
  auto make_Ru_3d = [](const cnt& m_cnt, const double shift_along_axis, const double z_shift, const double angle) {
    int n_atoms_in_cnt_unit_cell = m_cnt.pos_u_3d().n_rows;
    int total_number_of_atoms = m_cnt.pos_u_3d().n_rows * m_cnt.length_in_cnt_unit_cell();
    arma::mat all_atoms(total_number_of_atoms,3);

    for (int i=0; i<m_cnt.length_in_cnt_unit_cell(); i++)
    {
      all_atoms.rows(i*n_atoms_in_cnt_unit_cell,(i+1)*n_atoms_in_cnt_unit_cell-1) = i*m_cnt.pos_u_3d();
    }

    // make the cnt center at the middle
    double y_max = all_atoms.col(1).max();
    double y_min = all_atoms.col(1).min();
    all_atoms.col(1) -= ((y_max+y_min)/2.);

    // shift the center of the cnt axis along it's axis
    all_atoms.col(1) += shift_along_axis;

    // shift the atoms along the z axis
    all_atoms.col(2) +=  z_shift;

    // rotate by angle around the z axis
    for (unsigned int i=0; i<all_atoms.n_rows; i++)
    {
      double x = all_atoms(i,0)*std::cos(angle) - all_atoms(i,1)*std::sin(angle);
      double y = all_atoms(i,0)*std::sin(angle) + all_atoms(i,1)*std::cos(angle);
      all_atoms(i,0) = x;
      all_atoms(i,1) = y;
    }
    return all_atoms;
  };

    // make position of all atoms in the entire cnt length in 2d space of unrolled cnt
  auto make_Ru_2d = [](const cnt& m_cnt) {
    int n_atoms_in_cnt_unit_cell = m_cnt.pos_u_2d().n_rows;
    int total_number_of_atoms = m_cnt.pos_u_2d().n_rows * m_cnt.length_in_cnt_unit_cell();
    arma::mat all_atoms(total_number_of_atoms,2);

    for (int i=0; i<m_cnt.length_in_cnt_unit_cell(); i++)
    {
      all_atoms.rows(i*n_atoms_in_cnt_unit_cell,(i+1)*n_atoms_in_cnt_unit_cell-1) = i*m_cnt.pos_u_2d();
    }

    return all_atoms;
  };

  arma::mat i_Ru_3d = make_Ru_3d(*(pair.i.cnt_obj), shifts_along_axis[0], 0, 0);
  arma::mat f_Ru_3d = make_Ru_3d(*(pair.f.cnt_obj), shifts_along_axis[1], z_shift, angle);

  arma::mat i_Ru_2d = make_Ru_2d(*(pair.i.cnt_obj));
  arma::mat f_Ru_2d = make_Ru_2d(*(pair.f.cnt_obj));

  std::complex<double> J = 0;
  const std::complex<double> i1(0,1);

  // prebuild exponential factor for the inner loop
  arma::cx_vec f_exp(f_Ru_2d.n_rows);
  for (unsigned int j=0; j<f_Ru_2d.n_rows; j++)
  {
    f_exp(j) = std::exp(+i1*arma::dot(pair.f.ik_cm*pair.f.dk_l(),f_Ru_2d.row(j)));
  }

  for (unsigned int i=0; i<i_Ru_2d.n_rows; i++)
  {
    std::complex<double> i_exp = std::exp(-i1*arma::dot(pair.i.ik_cm*pair.i.dk_l(),i_Ru_2d.row(i)));    
    for (unsigned int j=0; j<f_Ru_2d.n_rows; j++)
    {
      J += i_exp*f_exp(j)/(arma::norm(i_Ru_3d.row(i)-f_Ru_3d.row(j)));
    }
  }
  return J;
}

// calculate first order transfer rate
double exciton_transfer::first_order(const double& z_shift, const std::array<double,2> axis_shifts, const double& theta, const bool& show_results)
{
  const cnt& donor = *_cnts[0];
  const cnt& acceptor = *_cnts[1];

  double min_energy = donor.A2_singlet().energy.min();

  const cnt::exciton_struct& d_exciton = donor.A2_singlet();
  const cnt::exciton_struct& a_exciton = acceptor.A2_singlet();

  // find lists of relevant states in donor and acceptor excitons
  std::vector<ex_state> d_relevant_states = get_relevant_states(d_exciton,min_energy);
  std::vector<ex_state> a_relevant_states = get_relevant_states(a_exciton,min_energy);

  double Z = 0;
  for (const auto& state:d_relevant_states)
  {
    Z += std::exp(-state.energy/(constants::kb*_temperature));
  }

  // match the states based on their energy
  std::vector<matching_states> state_pairs = match_states(d_relevant_states, a_relevant_states);

  double transfer_rate = 0;

  progress_bar prog(state_pairs.size(),"calculate first-order exciton transfer rate", not show_results);
  for (const auto& pair:state_pairs)
  { 
    prog.step();
    std::complex<double> Q = calculate_Q(pair);
    std::complex<double> J = calculate_J(pair, axis_shifts, z_shift, theta);
    double M = std::abs(Q*J)/std::sqrt(pair.i.cnt_obj->length_in_meter()*pair.f.cnt_obj->length_in_meter());
    transfer_rate += (2*constants::pi/constants::hb)*(std::exp(-pair.i.energy/(constants::kb*_temperature))/Z)*std::pow(M,2)*lorentzian(pair.i.energy-pair.f.energy);
  }

  if (show_results)
  {
    std::cout << "\n\n";
    std::cout << "cnt lengths: " << _cnts[0]->length_in_meter()*1.e9 << " [nm], " << _cnts[1]->length_in_meter()*1.e9 << " [nm]\n";
    std::cout << "center to center distance: " << z_shift*1.e9 << " [nm]\n";
    std::cout << "cnt1 radius: " << _cnts[0]->radius()*1.e9 << " [nm], cnt2 radius: " << _cnts[1]->radius()*1.e9 << " [nm]\n";
    std::cout << "wall to wall distance: " << (z_shift - _cnts[0]->radius() - _cnts[1]->radius())*1.e9 << " [nm]\n";
    std::cout << "theta: " << theta/constants::pi*180 << " [degrees]\n";
    std::cout << "axis shifts: " << axis_shifts[0]*1e9 << " [nm] and " << axis_shifts[1]*1e9 << " [nm]\n";
    std::cout << "exciton transfer rate: " << transfer_rate << "\n";
  }

  return transfer_rate;
}

// calculate first order transfer rate for varying angle
void exciton_transfer::calculate_first_order_vs_angle(const double& z_shift, const std::array<double,2> axis_shifts)
{
  int n_theta = 100;
  arma::vec theta_vec = arma::linspace<arma::vec>(0,constants::pi/2,n_theta);
  arma::vec transfer_rate(arma::size(theta_vec), arma::fill::zeros);

  progress_bar prog(n_theta, "first order transfer rate versus angle");
  for (int i=0; i<n_theta; i++)
  {
    prog.step();
    transfer_rate(i) = first_order(z_shift, axis_shifts, theta_vec(i));
  }

  // save the transfer rate
  std::string filename = _directory.path() / "first_order_transfer_rate_vs_angle.dat";
  transfer_rate.save(filename,arma::arma_ascii);

  // save theta vector
  filename = _directory.path() / "first_order_transfer_rate_vs_angle.theta.dat";
  theta_vec.save(filename,arma::arma_ascii);

  std::cout << "\n\n";
  std::cout << "cnt lengths: " << _cnts[0]->length_in_meter()*1.e9 << " [nm], " << _cnts[1]->length_in_meter()*1.e9 << " [nm]\n";
  std::cout << "center to center distance: " << z_shift*1.e9 << " [nm]\n";
  std::cout << "cnt1 radius: " << _cnts[0]->radius()*1.e9 << " [nm], cnt2 radius: " << _cnts[1]->radius()*1.e9 << " [nm]\n";
  std::cout << "wall to wall distance: " << (z_shift - _cnts[0]->radius() - _cnts[1]->radius())*1.e9 << " [nm]\n";
  std::cout << "axis shifts: " << axis_shifts[0]*1e9 << " [nm] and " << axis_shifts[1]*1e9 << " [nm]\n";
  std::cout << "max transfer rate: " << transfer_rate.max()/1e12 << " [1/ps] at " << theta_vec(transfer_rate.index_max())*180/constants::pi << " [degrees]\n";
  std::cout << "min transfer rate: " << transfer_rate.min()/1e12 << " [1/ps] at " << theta_vec(transfer_rate.index_min())*180/constants::pi << " [degrees]\n";
}

// calculate first order transfer rate for center to center distance
void exciton_transfer::calculate_first_order_vs_zshift(const std::array<double,2> axis_shifts, const double& theta)
{
  int n = 100;
  double z_shift_min = 1.5e-9;
  double z_shift_max = 1.5e-7;
  arma::vec z_shift = arma::linspace<arma::vec>(z_shift_min,z_shift_max,n);
  arma::vec transfer_rate(arma::size(z_shift), arma::fill::zeros);

  progress_bar prog(n, "first order transfer rate versus z_shift");
  for (int i=0; i<n; i++)
  {
    prog.step();
    transfer_rate(i) = first_order(z_shift(i), axis_shifts, theta);
  }

  // save the transfer rate
  std::string filename = _directory.path() / "first_order_transfer_rate_vs_zshift.dat";
  transfer_rate.save(filename,arma::arma_ascii);

  // save theta vector
  filename = _directory.path() / "first_order_transfer_rate_vs_zshift.distance.dat";
  z_shift.save(filename,arma::arma_ascii);

  std::cout << "\n\n";
  std::cout << "cnt lengths: " << _cnts[0]->length_in_meter()*1.e9 << " [nm], " << _cnts[1]->length_in_meter()*1.e9 << " [nm]\n";
  std::cout << "cnt1 radius: " << _cnts[0]->radius()*1.e9 << " [nm], cnt2 radius: " << _cnts[1]->radius()*1.e9 << " [nm]\n";
  std::cout << "theta: " << theta*180/constants::pi << " [degrees]\n";
  std::cout << "axis shifts: " << axis_shifts[0]*1e9 << " [nm] and " << axis_shifts[1]*1e9 << " [nm]\n";
  std::cout << "max transfer rate: " << transfer_rate.max() << " [1/s] at " << z_shift(transfer_rate.index_max())*1e9 << " [nm]\n";
  std::cout << "min transfer rate: " << transfer_rate.min() << " [1/s] at " << z_shift(transfer_rate.index_min())*1e9 << " [nm]\n";
}

// calculate first order transfer rate for varying axis shifts
void exciton_transfer::calculate_first_order_vs_axis_shift(const double z_shift, const double& theta)
{
  int n = 100;
  double axis_shift_min = 1.5e-9;
  double axis_shift_max = 1.5e-7;
  arma::vec shift_vec = arma::linspace<arma::vec>(axis_shift_min,axis_shift_max,n);
  arma::vec transfer_rate(arma::size(shift_vec), arma::fill::zeros);

  progress_bar prog(n, "first order transfer rate versus z_shift");
  for (int i=0; i<n; i++)
  {
    prog.step();
    std::array<double,2> axis_shifts = {0,shift_vec(i)};
    transfer_rate(i) = first_order(z_shift, axis_shifts, theta);
  }

  // save the transfer rate
  std::string filename = _directory.path() / "first_order_transfer_rate_vs_axis_shift.dat";
  transfer_rate.save(filename,arma::arma_ascii);

  // save theta vector
  filename = _directory.path() / "first_order_transfer_rate_vs_axis_shift.shift.dat";
  shift_vec.save(filename,arma::arma_ascii);

  std::cout << "\n\n";
  std::cout << "cnt lengths: " << _cnts[0]->length_in_meter()*1.e9 << " [nm], " << _cnts[1]->length_in_meter()*1.e9 << " [nm]\n";
  std::cout << "cnt1 radius: " << _cnts[0]->radius()*1.e9 << " [nm], cnt2 radius: " << _cnts[1]->radius()*1.e9 << " [nm]\n";
  std::cout << "theta: " << theta*180/constants::pi << " [degrees]\n";
  std::cout << "center to center distance: " << z_shift*1.e9 << " [nm]\n";
  std::cout << "max transfer rate: " << transfer_rate.max() << " [1/s] at " << shift_vec(transfer_rate.index_max())*1e9 << " [nm]\n";
  std::cout << "min transfer rate: " << transfer_rate.min() << " [1/s] at " << shift_vec(transfer_rate.index_min())*1e9 << " [nm]\n";
}