/**
cnt.cpp
Stores all relevant information for a carbon nanotube
*/

#include <iostream>
#include <numeric>
#include <armadillo>

#include "constants.h"
#include "cnt.h"

// set the output directory and the output file name
void cnt::process_command_line_args(int argc, char* argv[])
{
  namespace fs = std::experimental::filesystem;

  // first find the xml input file and open it
  fs::directory_entry xml_file;
  std::cout << "current path is " << fs::current_path() << std::endl;

  if (argc <= 1)
  {
    xml_file.assign("input.xml");
  }
  else
  {
    xml_file.assign(argv[1]);
  }

  if(fs::exists(xml_file))
  {
    std::cout << "input xml file found: " << xml_file.path() << std::endl;
  }
  else
  {
    std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
    std::exit(1);
  }

  if (!fs::is_regular_file(xml_file))
  {
    std::cout << "input xml file NOT found: " << xml_file.path() << std::endl;
    std::exit(1);
  }
  std::cout << std::endl;

  rapidxml::file<> xmlFile(xml_file.path().c_str()); //open file
  rapidxml::xml_document<> doc; //create xml object
  doc.parse<0>(xmlFile.data()); //parse contents of file
  rapidxml::xml_node<>* curr_node = doc.first_node(); //gets the node "Document" or the root nodes
  curr_node = curr_node->first_node();

  // get the sibling node with name sibling_name
  auto get_sibling = [](rapidxml::xml_node<>* node, const std::string& sibling_name)
  {
    if (node->name() == sibling_name)
    {
      return node;
    }
    auto next_node = node->next_sibling(sibling_name.c_str());
    if (next_node == 0)
    {
      next_node = node->previous_sibling(sibling_name.c_str());
      if (next_node == 0)
      {
        std::cout << "sibling not found: " << sibling_name.c_str() << std::endl;
        std::exit(1);
      }
    }
    return next_node;
  };

  // get name cnt name
  {
    curr_node = get_sibling(curr_node, "name");
    _name = trim_copy(curr_node->value());
    std::cout << "cnt name: '" << _name << "'\n";
  }

  // set the output_directory
  {
    curr_node = get_sibling(curr_node,"output_directory");

    std::string attr = curr_node->first_attribute("type")->value();
    // std::string path = trim_copy(curr_node->value());
    fs::path path = trim_copy(curr_node->value());
    if (attr == "absolute")
    {
      std::cout << "absolute directory format used!\n";
    }

    _directory.assign(path);
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
      std::exit(EXIT_FAILURE);
    }
  }

  // read chirality
  {
    curr_node = get_sibling(curr_node,"chirality");
    std::string chirality = curr_node->value();
    _n = std::stoi(chirality.substr(0,chirality.find(",")));
    _m = std::stoi(chirality.substr(chirality.find(",")+1));
    std::cout << "chirality: (" << _n << "," << _m << ")\n";
  }

  // length of cnt
  {
    curr_node = get_sibling(curr_node,"length");
    std::string units = curr_node->first_attribute("units")->value();
    if (units != "cnt_unit_cell")
    {
      std::cout << "cnt length is not in units '" << units << "'\n";
      std::exit(1);
    }
    _number_of_cnt_unit_cells = std::stoi(curr_node->value());
    std::cout << "length of cnt: " << _number_of_cnt_unit_cells << " unit cells.\n";
  }

  std::cout << std::endl;

};
// calculates position of atoms and reciprocal lattice vectors
void cnt::geometry()
{
  // unit vectors and reciprocal lattice vectors
  _a1 = arma::vec({_a_l*std::sqrt(3.0)/2.0, +_a_l/2.0});
  _a2 = arma::vec({_a_l*std::sqrt(3.0)/2.0, -_a_l/2.0});
  _b1 = arma::vec({1.0/sqrt(3.0)*2.0*constants::pi/_a_l, +2.0*constants::pi/_a_l});
  _b2 = arma::vec({1.0/sqrt(3.0)*2.0*constants::pi/_a_l, -2.0*constants::pi/_a_l});
  _b2 = arma::vec({_a_l*std::sqrt(3.0)/2.0, _a_l/2.0});

  _a1.print("a1:");
  _a2.print("a2:");
  _b1.print("b1:");
  _b2.print("b2:");


	_aCC_vec = 1.0/3.0*(_a1+_a2);

  _aCC_vec.print("aCC vector:");

	// calculate chirality and translational vectors of CNT unit cell
	_ch_vec = double(_n) * _a1 + double(_m) * _a2;
	_ch_len = arma::norm(_ch_vec,2);

  _ch_vec.print("chirality vector:");
  std::cout << "ch_vec length:\n   " << _ch_len << std::endl;


	_radius = _ch_len/2.0/constants::pi;

	int dR = std::gcd(2*_n+_m,_n+2*_m);
	int t1 = +(2*_m+_n)/dR;
	int t2 = -(2*_n+_m)/dR;

	_t_vec = double(t1)*_a1 + double(t2)*_a2;
  _t_vec.print("t_vec:");

	_Nu = 2*(std::pow(_n,2)+std::pow(_m,2)+_n*_m)/dR;

	// rotate basis vectors so that ch_vec is along the x_axis
	double cos_theta = _ch_vec(0)/arma::norm(_ch_vec);
	double sin_theta = _ch_vec(1)/arma::norm(_ch_vec);
	arma::mat rot = {{+cos_theta, +sin_theta},
                   {-sin_theta, +cos_theta}}; // rotation matrix

	_ch_vec = rot*_ch_vec;
	_t_vec = rot*_t_vec;
	_a1 = rot*_a1;
	_a2 = rot*_a2;
	_b1 = rot*_b1;
	_b2 = rot*_b2;
	_aCC_vec = rot*_aCC_vec;

  _ch_vec.print("chirality vector:");
  _t_vec.print("t_vec:");


	// calculate reciprocal lattice of CNT
	_K1 = (-double(t2)*_b1 + double(t1)*_b2)/(double(_Nu));
	_K2 = (double(_m)*_b1-double(_n)*_b2)/(double(_Nu));
	_dk_l = _K2/(double(_number_of_cnt_unit_cells));
  _nk = _number_of_cnt_unit_cells;

	// calculate positions of atoms in the cnt unit cell
	_pos_a = arma::mat(_Nu,2,arma::fill::zeros);
	_pos_b = arma::mat(_Nu,2,arma::fill::zeros);

	int k = 0;

	for (int i=0; i<=t1+_n; i++)
	{
		for (int j=t2; j<=_m; j++)
		{
			bool flag1 = double(t2*i)/(double)t1 <= double(j);
			bool flag2 = double(_m*i)/(double)_n >= double(j);
			bool flag3 = double(t2*(i-_n))/double(t1) > double(j-_m);
			bool flag4 = double(_m*(i-t1))/double(_n) < double(j-t2);

			if(flag1 && flag2 && flag3 && flag4)
			{
        _pos_a.row(k) = double(i)*_a1.t() + double(j)*_a2.t();
        _pos_b.row(k) = _pos_a.row(k) + _aCC_vec.t();

				if(_pos_a(k,0) > _ch_vec(0))
          _pos_a(k,0) -= _ch_vec(0);
				if(_pos_a(k,0) < 0.0)
          _pos_a(k,0) += _ch_vec(0);
				if(_pos_a(k,1) > _ch_vec(1))
          _pos_a(k,1) -= _ch_vec(1);
				if(_pos_a(k,1) < 0.0)
          _pos_a(k,1) += _ch_vec(1);

				if(_pos_b(k,0) > _ch_vec(0))
          _pos_b(k,0) -= _ch_vec(0);
				if(_pos_b(k,0) < 0.0)
          _pos_b(k,0) += _ch_vec(0);
				if(_pos_b(k,1) > _ch_vec(1))
          _pos_b(k,1) -= _ch_vec(1);
				if(_pos_b(k,1) < 0.0)
          _pos_b(k,1) += _ch_vec(1);

				k++;
			}
		}
	}

  _pos_a.print("pos_a:");
  _pos_b.print("pos_b:");

	if (k != _Nu)
	{
		std::cout << "error in finding position of atoms in cnt unit cell!!!" << std::endl;
		std::cout << "Nu = " << _Nu << "  ,  k = " << k << std::endl;
		exit(1);
	}

	// calculate distances between atoms in a warped cnt unit cell.
	_pos_aa = arma::mat(_Nu,2,arma::fill::zeros);
	_pos_ab = arma::mat(_Nu,2,arma::fill::zeros);
	_pos_ba = arma::mat(_Nu,2,arma::fill::zeros);
	_pos_bb = arma::mat(_Nu,2,arma::fill::zeros);

	for (int i=0; i<_Nu; i++)
	{
    _pos_aa.row(i) = _pos_a.row(i)-_pos_a.row(0);
    _pos_ab.row(i) = _pos_a.row(i)-_pos_b.row(0);
    _pos_ba.row(i) = _pos_b.row(i)-_pos_a.row(0);
    _pos_bb.row(i) = _pos_b.row(i)-_pos_b.row(0);

		if(_pos_aa(i,0) > _ch_vec(0)/2)
      _pos_aa(i,0) -= _ch_vec(0);
		if(_pos_ab(i,0) > _ch_vec(0)/2)
      _pos_ab(i,0) -= _ch_vec(0);
		if(_pos_ba(i,0) > _ch_vec(0)/2)
      _pos_ba(i,0) -= _ch_vec(0);
		if(_pos_bb(i,0) > _ch_vec(0)/2)
      _pos_bb(i,0) -= _ch_vec(0);
	}

	// put position of all atoms in a single variable in 2d space(unrolled graphene sheet)
	_pos_2d = arma::mat(2*_Nu,2,arma::fill::zeros);
  _pos_2d(arma::span(0,_Nu-1),arma::span::all) = _pos_a;
  _pos_2d(arma::span(_Nu,2*_Nu-1),arma::span::all) = _pos_b;

	// calculate position of all atoms in the 3d space (rolled graphene sheet)
	_pos_3d = arma::mat(2*_Nu,3,arma::fill::zeros);
	for (int i=0; i<_pos_3d.n_rows; i++)
	{
		_pos_3d(i,0) = _radius*cos(_pos_2d(i,0)/_radius);
		_pos_3d(i,1) = _pos_2d(i,1);
		_pos_3d(i,2) = _radius*sin(_pos_2d(i,0)/_radius);
	}

	//make 3d t_vec
	_t_vec_3d = arma::vec(3,arma::fill::zeros);
	_t_vec_3d(1) = _t_vec(1);


	// save coordinates of atoms in 2d space
  std::string filename = _directory.path().string() + _name + ".pos_2d.dat";
  _pos_2d.save(filename, arma::arma_ascii);

  // save coordinates of atoms in 3d space
  filename = _directory.path().string() + _name + ".pos_3d.dat";
  _pos_3d.save(filename, arma::arma_ascii);

}
// calculate electron dispersion energies using full unit cell (2*Nu atoms)
void cnt::electron_full()
{

	// make the list of 1st nearest neighbor atoms
	arma::umat nn_list(2*_Nu,3,arma::fill::zeros); // contains index of the nearest neighbor atom
	arma::imat nn_tvec_index(2*_Nu,3,arma::fill::zeros); // contains the index of the cnt unit cell that the nearest neigbor atom is in.
	for (int i=0; i<_pos_3d.n_rows; i++)
	{
		int k=0;
		for (int j=0; j<_pos_3d.n_rows; j++)
		{
			for (int l=-1; l<=1; l++)
			{
        double dR = arma::norm(_pos_3d.row(i)-_pos_3d.row(j)+double(l)*_t_vec_3d.t());
				if ( (i!=j) && (dR<(1.4*_a_cc)) )
				{
					nn_list(i,k) = j;
					nn_tvec_index(i,k) = l;
					k++;
				}
			}
		}
    if (k != 3)
    {
      std::cout << "error: nearest neighbors partially found!!!\n";
      std::exit(1);
    }
	}

	arma::cx_mat H(2*_Nu, 2*_Nu, arma::fill::zeros);
	arma::cx_mat S(2*_Nu, 2*_Nu, arma::fill::zeros);
  arma::vec E;
  arma::cx_mat C;

	int NK = _nk;

	_el_energy_full = arma::mat(2*_Nu, NK, arma::fill::zeros);
	_el_psi_full = arma::cx_cube(2*_Nu, 2*_Nu, NK, arma::fill::zeros);

	double t_len = arma::norm(_t_vec_3d);

	for (int n=0; n<NK; n++)
	{
		double wave_vec = double(n-_nk/2)*arma::norm(_dk_l);

		H.zeros();
		S.zeros();

		for (int i=0; i<2*_Nu; i++)
		{
			H(i,i) = (_e2p,0.e0);
			for (int k=0; k<3; k++)
			{
				int j = nn_list(i,k);
				int l = nn_tvec_index(i,k);

				H(i,j) += arma::cx_double(_t0,0.e0)*exp(arma::cx_double(0.0,wave_vec*double(l)*t_len));
				S(i,j) += arma::cx_double(_s0,0.e0)*exp(arma::cx_double(0.0,wave_vec*double(l)*t_len));
			}
		}

		arma::eig_sym(E, C, H);

		// fix the phase of the eigen vectors
		for (int i=0; i<C.n_cols; i++)
		{
			arma::cx_double phi = std::conj(C(0,i))/std::abs(C(0,i));
			for (int j=0; j<C.n_rows; j++)
			{
				C(j,i) *= phi;
			}
		}

    _el_energy_full.col(n) = E;
    _el_psi_full.slice(n) = C;

	}

  // save electron energy bands using full Brillouine zone
  std::string filename = _directory.path().string() + _name + ".el_energy_full.dat";
  _el_energy_full.save(filename, arma::arma_ascii);

  // // save electron wavefunctions using full Brillouine zone
  // filename = _directory.path().string() + _name + ".el_psi_full.dat";
  // _el_psi_full.save(filename, arma::arma_ascii);

}



//
//
// void cnt::dielectric()
// {
//
//
// 	// calculate v_q
// 	int Nq = nk;
// 	nr::mat3d_complex v_q(2*Nu, 2*Nu, Nq, nr::cmplx(0.0));
// 	nr::vec_doub wave_vec(Nq, 0.0);
//
// 	double t_len = t_vec_3d.norm2();
// 	double factor = nr::sqr(4*constants::pi*constants::eps0*Upp/nr::sqr(constants::q0));
// 	for(int iq = 0; iq<Nq; iq++)
// 	{
// 		wave_vec(iq) = double(iq-Nq/2)*dk_l.norm2();
// 		for (int i=-int(number_of_cnt_unit_cells/2); i<int(number_of_cnt_unit_cells/2); i++)
// 		{
// 			for (int b=0; b<2*Nu; b++)
// 			{
// 				for (int bp=0; bp<2*Nu; bp++)
// 				{
// 					double dx = pos_3d(b,0) - (pos_3d(bp,0) + double(i)*t_vec_3d(0));
// 					double dy = pos_3d(b,1) - (pos_3d(bp,1) + double(i)*t_vec_3d(1));
// 					double dz = pos_3d(b,2) - (pos_3d(bp,2) + double(i)*t_vec_3d(2));
// 					double dR2 = dx*dx+dy*dy+dz*dz;
//
// 					double I = Upp/sqrt(1+factor*dR2);
// 					v_q(b,bp,iq) += nr::cmplx(1.0/double(number_of_cnt_unit_cells))*exp(nr::cmplx(0.0,wave_vec(iq)*double(i)*t_len))*I;
// 				}
// 			}
// 		}
// 	}
//
// 	// open file to save v_q
// 	std::ofstream file_vq_real, file_vq_imag;
// 	file_vq_real.open(name+".vq.real.dat", std::ios::app);
// 	file_vq_imag.open(name+".vq.imag.dat", std::ios::app);
//
// 	file_vq_real << std::scientific;
// 	file_vq_real << std::showpos;
// 	file_vq_imag << std::scientific;
// 	file_vq_imag << std::showpos;
// 	for (int iq=0; iq<v_q.dim3(); iq++)
// 	{
// 		file_vq_real << wave_vec(iq) << "\t";
// 		file_vq_imag << wave_vec(iq) << "\t";
// 	}
// 	file_vq_real << "\n";
// 	file_vq_imag << "\n";
// 	for (int b=0; b<v_q.dim1(); b++)
// 	{
// 		for (int bp=0; bp<v_q.dim2(); bp++)
// 		{
// 			for (int iq=0; iq<v_q.dim3(); iq++)
// 			{
// 				file_vq_real << std::real(v_q(b,bp,iq)/constants::eV) << "\t";
// 				file_vq_imag << std::imag(v_q(b,bp,iq)/constants::eV) << "\t";
// 			}
// 			file_vq_real << "\n";
// 			file_vq_imag << "\n";
// 		}
// 	}
// 	file_vq_real.close();
// 	file_vq_imag.close();
//
//
// 	// // calculate polarization (PI)
// 	// nr::vec_doub PI(Nq,0.0);
// 	// for (int i=0; i<Nq; i++)
// 	// {
// 	// 	int iq = i-Nq/2;
// 	// 	nr::cmplx num1, num2;
// 	// 	double denom1, denom2;
//
// 	// 	for (int ik=0; ik<nk; ik++)
// 	// 	{
// 	// 		int ik_p = iq+ik;
// 	// 		while(ik_p>=nk)	ik_p = ik_p-nk;
// 	// 		while(ik_p<0) ik_p = ik_p+nk;
//
// 	// 		for (int alpha_v=0; alpha_v<Nu; alpha_v++)
// 	// 		{
// 	// 			for (int alpha_c=Nu; alpha_c<2*Nu; alpha_c++)
// 	// 			{
// 	// 				num1 = nr::cmplx(0.0,0.0);
// 	// 				num2 = nr::cmplx(0.0,0.0);
// 	// 				for(int b=0; b<2*Nu; b++)
// 	// 				{
// 	// 					num1 += std::conj(el_psi(b,alpha_v,ik)) * el_psi(b,alpha_c,ik_p);
// 	// 					num2 += std::conj(el_psi(b,alpha_c,ik)) * el_psi(b,alpha_v,ik_p);
// 	// 				}
// 	// 				denom1 = el_energy(alpha_c,ik_p)-el_energy(alpha_v,ik);
// 	// 				denom2 = el_energy(alpha_c,ik)-el_energy(alpha_v,ik_p);
//
// 	// 				PI(i) += std::norm(num1)/denom1 + std::norm(num2)/denom2;
// 	// 			}
// 	// 		}
// 	// 	}
//
// 	// 	std::cout << "calculating dielectric function: PI(" << i << ") = " << PI(i) << std::endl;
// 	// }
//
// 	// PI = 2.0*PI;
//
// 	// // open file to save polarization (PI)
// 	// std::ofstream file_pi;
// 	// file_pi.open(name+".pi.dat", std::ios::app);
//
// 	// file_pi << std::scientific;
// 	// file_pi << std::showpos;
// 	// for (int i=0; i<wave_vec.size(); i++)
// 	// {
// 	// 	file_pi << wave_vec(i) << "\t";
// 	// }
// 	// file_pi << "\n";
// 	// for (int i=0; i<PI.size(); i++)
// 	// {
// 	// 	file_pi << PI(i) << "\t";
// 	// }
// 	// file_pi.close();
//
// 	// // calculate dielectric function (epsilon)
// 	// epsilon.assign(Nq,nr::cmplx(0));
// 	// for (int i=0; i<epsilon.size(); i++)
// 	// {
// 	// 	nr::cmplx sum = nr::cmplx(0);
// 	// 	for (int b=0; b<2*Nu; b++)
// 	// 	{
// 	// 		for (int bp=0; bp<2*Nu; bp++)
// 	// 		{
// 	// 			sum += v_q(b,bp,i);
// 	// 		}
// 	// 	}
// 	// 	sum = sum/nr::cmplx((2*Nu)*(2*Nu));
// 	// 	epsilon(i) = nr::cmplx(1)+sum*nr::cmplx(PI(i));
// 	// }
//
// 	// // open file to save epsilon
// 	// std::ofstream file_eps;
// 	// file_eps.open(name+".epsilon.dat", std::ios::app);
//
// 	// file_eps << std::scientific;
// 	// file_eps << std::showpos;
// 	// for (int i=0; i<wave_vec.size(); i++)
// 	// {
// 	// 	file_eps << wave_vec(i) << "\t";
// 	// }
// 	// file_eps << "\n";
// 	// for (int i=0; i<epsilon.size(); i++)
// 	// {
// 	// 	file_eps << std::real(epsilon(i)) << "\t";
// 	// }
// 	// file_eps << "\n";
// 	// for (int i=0; i<epsilon.size(); i++)
// 	// {
// 	// 	file_eps << std::imag(epsilon(i)) << "\t";
// 	// }
// 	// file_eps.close();
// }
//
//
// void cnt::coulomb_int()
// {
//
//
// 	// calculate v_q
// 	const int Nq = nk;
// 	nr::mat3d_complex v_q(2*Nu, 2*Nu, Nq, nr::cmplx(0.0));
// 	nr::vec_doub wave_vec(Nq, 0.0);
//
// 	double t_len = t_vec_3d.norm2();
// 	double factor = nr::sqr(4*constants::pi*constants::eps0*Upp/nr::sqr(constants::q0));
// 	for(int iq = 0; iq<Nq; iq++)
// 	{
// 		wave_vec(iq) = double(iq)*dk_l.norm2();
// 		for (int i=-int(number_of_cnt_unit_cells/2); i<int(number_of_cnt_unit_cells/2); i++)
// 		{
// 			for (int b=0; b<2*Nu; b++)
// 			{
// 				for (int bp=0; bp<2*Nu; bp++)
// 				{
// 					double dx = pos_3d(b,0) - (pos_3d(bp,0) + double(i)*t_vec_3d(0));
// 					double dy = pos_3d(b,1) - (pos_3d(bp,1) + double(i)*t_vec_3d(1));
// 					double dz = pos_3d(b,2) - (pos_3d(bp,2) + double(i)*t_vec_3d(2));
// 					double dR2 = dx*dx+dy*dy+dz*dz;
//
// 					double I = Upp/sqrt(1+factor*dR2);
// 					v_q(b,bp,iq) += nr::cmplx(1.0/double(number_of_cnt_unit_cells))*exp(nr::cmplx(0.0,wave_vec(iq)*double(i)*t_len))*I;
// 				}
// 			}
// 		}
// 	}
//
// 	const int Nb = 4;
//
// 	// direct interaction matrix
// 	nr::mat_complex V_dir(Nb*Nb*nk,Nb*Nb*nk,nr::cmplx(0));
// 	nr::mat_complex V_xch(Nb*Nb*nk,Nb*Nb*nk,nr::cmplx(0));
// 	nr::mat_complex dE(Nb*Nb*nk,Nb*Nb*nk,nr::cmplx(0));
//
// 	nr::mat_complex tmp_dir(2*Nu,2*Nu,nr::cmplx(0));
// 	nr::mat_complex tmp_xch(2*Nu,2*Nu,nr::cmplx(0));
//
// 	const int iKcm=0; // center of mass wave vector
//
// 	for (int ic=0; ic<Nb; ic++)
// 	{
// 		int ac = Nu+ic;
// 		std::cout << "ic = " << ic << "\n";
// 		for (int icp=0; icp<Nb; icp++)
// 		{
// 			int acp = Nu+icp;
// 			for (int iv=0; iv<Nb; iv++)
// 			{
// 				int av = Nu-1-iv;
// 				for (int ivp=0; ivp<Nb; ivp++)
// 				{
// 					int avp = Nu-1-ivp;
// 					for (int ikc=0; ikc<nk; ikc++)
// 					{
// 						int ikv = int(ikc+iKcm);
// 						while(ikv>=nk)	ikv -= nk;
// 						while(ikv<0) ikv += nk;
//
// 						// int iq_xch = int(ikc-ikv);
// 						// while(iq_xch>=Nq) iq_xch -= Nq;
// 						// while(iq_xch<0) iq_xch += Nq;
//
// 						// for (int b=0; b<2*Nu; b++)
// 						// {
// 						// 	for (int bp=0; bp<2*Nu; bp++)
// 						// 	{
// 						// 		tmp_dir(b,bp) = std::conj(el_psi(b,ac,ikc))*el_psi(bp,av,ikv)/kappa;
// 						// 		tmp_xch(b,bp) = std::conj(el_psi(b,ac,ikc))*el_psi(b,av,ikv)*v_q(b,bp,iq_xch);
// 						// 	}
// 						// }
//
// 						// int row = (nk)*((Nb)*ic+iv)+ikc;
//
// 						dE((nk)*((Nb)*ic+iv)+ikc,(nk)*((Nb)*ic+iv)+ikc) = nr::cmplx(el_energy(ac,ikc)-el_energy(av,ikv));
//
// 						for (int ikcp=0; ikcp<nk; ikcp++)
// 						{
// 							int ikvp = int(ikcp+iKcm);
// 							while(ikv>=nk)	ikv -= nk;
// 							while(ikv<0) ikv += nk;
//
// 							int iq_dir = ikc-ikcp;
// 							while(iq_dir>=Nq) iq_dir = iq_dir-Nq;
// 							while(iq_dir<0) iq_dir = iq_dir+Nq;
//
// 							int iq_xch = int(ikc-ikv);
// 							while(iq_xch>=Nq) iq_xch -= Nq;
// 							while(iq_xch<0) iq_xch += Nq;
//
// 							// int col = (nk)*((Nb)*icp+ivp)+ikcp;
//
// 							for (int b=0; b<2*Nu; b++)
// 							{
// 								for (int bp=0; bp<2*Nu; bp++)
// 								{
// 									V_dir((nk)*((Nb)*ic+iv)+ikc,(nk)*((Nb)*icp+ivp)+ikcp) += std::conj(el_psi(b,ac,ikc))*el_psi(b,acp,ikcp)*el_psi(bp,av,ikv)*std::conj(el_psi(bp,avp,ikvp))*v_q(b,bp,iq_dir)/kappa;
// 									V_xch((nk)*((Nb)*ic+iv)+ikc,(nk)*((Nb)*icp+ivp)+ikcp) += std::conj(el_psi(b,ac,ikc))*el_psi(b,av,ikv)*el_psi(bp,acp,ikcp)*std::conj(el_psi(bp,avp,ikvp))*v_q(b,bp,iq_xch);
//
// 									// V_dir(row,col) += tmp_dir(b,bp)*el_psi(b,acp,ikcp)*std::conj(el_psi(bp,avp,ikvp))*v_q(b,bp,iq_dir);
// 									// V_xch(row,col) += tmp_xch(b,bp)*el_psi(bp,acp,ikcp)*std::conj(el_psi(bp,avp,ikvp));
//
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
//
//
// 	// open file to save V_dir
// 	std::ofstream file_V_dir_real, file_V_dir_imag;
// 	file_V_dir_real.open(name+".v_dir.real.dat", std::ios::app);
// 	file_V_dir_imag.open(name+".v_dir.imag.dat", std::ios::app);
//
// 	file_V_dir_real << std::scientific << std::showpos;
// 	file_V_dir_imag << std::scientific << std::showpos;
//
// 	file_V_dir_real << double(0.0) << "\t";
// 	file_V_dir_imag << double(0.0) << "\t";
// 	for (int ikcp=0; ikcp<V_dir.dim2(); ikcp++)
// 	{
// 		double wave_vec = double(ikcp)*dk_l.norm2();
// 		file_V_dir_real << wave_vec << "\t";
// 		file_V_dir_imag << wave_vec << "\t";
// 	}
// 	file_V_dir_real << "\n";
// 	file_V_dir_imag << "\n";
//
// 	for (int ikc=0; ikc<V_dir.dim1(); ikc++)
// 	{
// 		double wave_vec = double(ikc)*dk_l.norm2();
// 		file_V_dir_real << wave_vec << "\t";
// 		file_V_dir_imag << wave_vec << "\t";
//
// 		for (int ikcp=0; ikcp<V_dir.dim2(); ikcp++)
// 		{
// 			file_V_dir_real << std::real(V_dir(ikc,ikcp)/constants::eV) << "\t";
// 			file_V_dir_imag << std::imag(V_dir(ikc,ikcp)/constants::eV) << "\t";
// 		}
// 		file_V_dir_real << "\n";
// 		file_V_dir_imag << "\n";
// 	}
// 	file_V_dir_real.close();
// 	file_V_dir_imag.close();
//
// 	// open file to save V_xch
// 	std::ofstream file_V_xch_real, file_V_xch_imag;
// 	file_V_xch_real.open(name+".v_xch.real.dat", std::ios::app);
// 	file_V_xch_imag.open(name+".v_xch.imag.dat", std::ios::app);
//
// 	file_V_xch_real << std::scientific << std::showpos;
// 	file_V_xch_imag << std::scientific << std::showpos;
//
// 	file_V_xch_real << double(0.0) << "\t";
// 	file_V_xch_imag << double(0.0) << "\t";
// 	for (int ikcp=0; ikcp<V_xch.dim2(); ikcp++)
// 	{
// 		double wave_vec = double(ikcp)*dk_l.norm2();
// 		file_V_xch_real << wave_vec << "\t";
// 		file_V_xch_imag << wave_vec << "\t";
// 	}
// 	file_V_xch_real << "\n";
// 	file_V_xch_imag << "\n";
//
// 	for (int ikc=0; ikc<V_xch.dim1(); ikc++)
// 	{
// 		double wave_vec = double(ikc)*dk_l.norm2();
// 		file_V_xch_real << wave_vec << "\t";
// 		file_V_xch_imag << wave_vec << "\t";
//
// 		for (int ikcp=0; ikcp<V_xch.dim2(); ikcp++)
// 		{
// 			file_V_xch_real << std::real(V_xch(ikc,ikcp)/constants::eV) << "\t";
// 			file_V_xch_imag << std::imag(V_xch(ikc,ikcp)/constants::eV) << "\t";
// 		}
// 		file_V_xch_real << "\n";
// 		file_V_xch_imag << "\n";
// 	}
// 	file_V_xch_real.close();
// 	file_V_xch_imag.close();
//
//
// 	// solve for exciton energy
// 	nr::mat_complex kernel(dE-V_dir);
// 	// nr::mat_complex kernel(dE+nr::cmplx(2)*V_xch-V_dir);
//
// 	nr::mat_complex ex_psi(Nb*Nb*nk,Nb*Nb*nk);
// 	nr::vec_doub ex_energy(Nb*Nb*nk);
// 	nr::eig_sym(ex_energy, ex_psi, kernel);
//
// 	// nr::mat_complex ex_psi(Nb*Nb*nk,Nb*Nb*nk,nr::cmplx(0));
// 	// nr::vec_doub ex_energy(Nb*Nb*nk,0.0);
// 	// nr::eig_sym_selected(ex_energy, ex_psi, kernel, 20);
//
//
// 	// open file to save exciton energies
// 	std::ofstream file_ex_energy;
// 	file_ex_energy.open(name+".ex_energy.dat", std::ios::app);
// 	file_ex_energy << std::scientific << std::showpos;
// 	for (int i=0; i<ex_energy.size(); i++)
// 	{
// 		file_ex_energy << ex_energy(i)/constants::eV << "\t" << std::real(dE(i,i))/constants::eV << "\n";
// 	}
// 	file_ex_energy.close();
//
//
// 	// open file to save exciton energies
// 	std::ofstream file_ex_psi;
// 	file_ex_psi.open(name+".ex_psi.dat", std::ios::app);
// 	file_ex_psi << std::scientific << std::showpos;
// 	for (int i=0; i<ex_psi.dim1(); i++)
// 	{
// 		for (int j=0; j<ex_psi.dim2(); j++)
// 		{
// 			file_ex_psi << std::real(ex_psi(i,j)) << "\t";
// 		}
// 		file_ex_psi << "\n";
// 	}
// 	file_ex_psi.close();
//
// }
