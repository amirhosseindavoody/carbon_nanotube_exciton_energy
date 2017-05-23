/**
cnt.cpp
Stores all relevant information for a carbon nanotube
*/

#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>	// for opening the output files

#include "cnt.h"
#include "constants.h"
#include "write_log.h"
#include "nr.h"

cnt::cnt(const std::string &in_name, const int in_n, const int in_m, const int in_length)
{
	name = in_name;
	n = in_n;
	m = in_m;
	number_of_cnt_unit_cells = in_length;
	nk = number_of_cnt_unit_cells;

	a1.assign(2,0.0);
	a2.assign(2,0.0);
	b1.assign(2,0.0);
	b2.assign(2,0.0);

}

void cnt::geometry()
{

	// unit vectors and reciprocal lattice vectors
	a1[0] = a_l*sqrt(3.0e+0)/2.0e+0;
	a1[1] = a_l/2.0e+0;

	a2[0] = a_l*sqrt(3.0e+0)/2.0e+0;
	a2[1] = -a_l/2.0e+0;
	
	b1[0] = 1.0/sqrt(3.0)*2.0*constants::pi/a_l;
	b1[1] = +2.0*constants::pi/a_l;
	
	b2[0] = 1.0/sqrt(3.0)*2.0*constants::pi/a_l;
	b2[1] = -2.0*constants::pi/a_l;

	aCC_vec = 1.0/3.0*(a1+a2);

	// calculate chirality and translational vectors of CNT unit cell
	ch_vec = (double)n * a1 + (double)m * a2;
	ch_len = a_l*sqrt(pow((double)n,2)+pow((double)m,2)+(double)n*m);

	radius = ch_len/2.0/constants::pi;

	int dR = nr::gcd(2*n+m,n+2*m);
	int t1 = +(2*m+n)/dR;
	int t2 = -(2*n+m)/dR;
	t_vec = (double)t1*a1 + (double)t2*a2;

	Nu = 2*(pow(n,2)+pow(m,2)+n*m)/dR;

	// rotate basis vectors so that ch_vec is along the x_axis
	double cos_theta = ch_vec[0]/ch_vec.norm2();
	double sin_theta = ch_vec[1]/ch_vec.norm2();
	nr::mat_doub rot(2,2); // rotation matrix
	rot(0,0) = cos_theta;
	rot(1,0) = -sin_theta;
	rot(0,1) = sin_theta;
	rot(1,1) = cos_theta;

	ch_vec = rot*ch_vec;
	t_vec = rot*t_vec;
	a1 = rot*a1;
	a2 = rot*a2;
	b1 = rot*b1;
	b2 = rot*b2;
	aCC_vec = rot*aCC_vec;
	
	// calculate reciprocal lattice of CNT
	K1 = (-(double)t2*b1 + (double)t1*b2)/((double)Nu);
	K2 = ((double)m*b1-(double)n*b2)/((double)Nu);
	dk_l = K2/((double)number_of_cnt_unit_cells);

	// calculate positions of atoms in the cnt unit cell
	pos_a.assign(Nu,2,0.0);
	pos_b.assign(Nu,2,0.0);

	int k =0;

	for (int i=0; i<=t1+n; i++)
	{
		for (int j=t2; j<=m; j++)
		{
			bool flag1 = (double)(t2*i)/(double)t1 <= (double)j;
			bool flag2 = (double)(m*i)/(double)n >= (double)j;
			bool flag3 = (double)(t2*(i-n))/(double)(t1) > (double)(j-m);
			bool flag4 = (double)(m*(i-t1))/(double)(n) < (double)(j-t2);

			if(flag1 && flag2 && flag3 && flag4)
			{
				pos_a(k,0) = (double)i*a1(0) + (double)j*a2(0);
				pos_a(k,1) = (double)i*a1(1) + (double)j*a2(1);
				pos_b(k,0) = pos_a(k,0)+aCC_vec(0);
				pos_b(k,1) = pos_a(k,1)+aCC_vec(1);

				if(pos_a(k,0) > ch_vec(0))	pos_a(k,0) -= ch_vec(0);
				if(pos_a(k,0) < 0.0)	pos_a(k,0) += ch_vec(0);
				if(pos_a(k,1) > ch_vec(1))	pos_a(k,1) -= ch_vec(1);
				if(pos_a(k,1) < 0.0)	pos_a(k,1) += ch_vec(1);

				if(pos_b(k,0) > ch_vec(0))	pos_b(k,0) -= ch_vec(0);
				if(pos_b(k,0) < 0.0)	pos_b(k,0) += ch_vec(0);
				if(pos_b(k,1) > ch_vec(1))	pos_b(k,1) -= ch_vec(1);
				if(pos_b(k,1) < 0.0)	pos_b(k,1) += ch_vec(1);
				
				k++;
			}
		}
	}

	if (k != Nu) 
	{
		std::cout << "error in finding position of atoms in cnt unit cell!!!" << std::endl;
		std::cout << "Nu = " << Nu << "  ,  k = " << k << std::endl;
		exit(1);
	}

	// // calculate distances between atoms in a warped cnt unit cell.
	// pos_aa.assign(Nu,2,0.0);
	// pos_ab.assign(Nu,2,0.0);
	// pos_ba.assign(Nu,2,0.0);
	// pos_bb.assign(Nu,2,0.0);

	// for (int i=0; i<Nu; i++)
	// {
	// 	pos_aa(i,0) = pos_a(i,0)-pos_a(0,0);
	// 	pos_aa(i,1) = pos_a(i,1)-pos_a(0,1);

	// 	pos_ab(i,0) = pos_a(i,0)-pos_b(0,0);
	// 	pos_ab(i,1) = pos_a(i,1)-pos_b(0,1);

	// 	pos_ba(i,0) = pos_b(i,0)-pos_a(0,0);
	// 	pos_ba(i,1) = pos_b(i,1)-pos_a(0,1);

	// 	pos_bb(i,0) = pos_b(i,0)-pos_b(0,0);
	// 	pos_bb(i,1) = pos_b(i,1)-pos_b(0,1);

	// 	if(pos_aa(i,0) > ch_vec(0)/2)	pos_aa(i,0) -= ch_vec(0);
	// 	if(pos_ab(i,0) > ch_vec(0)/2)	pos_ab(i,0) -= ch_vec(0);
	// 	if(pos_ba(i,0) > ch_vec(0)/2)	pos_ba(i,0) -= ch_vec(0);
	// 	if(pos_bb(i,0) > ch_vec(0)/2)	pos_bb(i,0) -= ch_vec(0);
	// }

	// put position of all atoms in a single variable in 2d space(unrolled graphene sheet)
	pos_2d.assign(2*Nu,2,0.0);
	for (int i=0; i<Nu; i++)
	{
		pos_2d(i,0) = pos_a(i,0);
		pos_2d(i+Nu,0) = pos_b(i,0);
		pos_2d(i,1) = pos_a(i,1);
		pos_2d(i+Nu,1) = pos_b(i,1);
	}

	// calculate position of all atoms in the 3d space (rolled graphene sheet)
	pos_3d.assign(2*Nu,3,0.0);
	for (int i=0; i<pos_3d.dim1(); i++)
	{
		pos_3d(i,0) = radius*cos(pos_2d(i,0)/radius);
		pos_3d(i,1) = pos_2d(i,1);
		pos_3d(i,2) = radius*sin(pos_2d(i,0)/radius);
	}

	//make 3d t_vec
	t_vec_3d.assign(3,0.0);
	t_vec_3d(1) = t_vec(1);


	// save coordinates of atoms in 2d space
	std::ofstream file;
	file.open(name+".pos_2d.dat", std::ios::trunc);

	if (file.is_open())
	{
		file << std::scientific;
		file << std::showpos;
		file << pos_2d.dim1() << "\t" << pos_2d.dim2() << std::endl;
		for (int i=0; i<pos_2d.dim1(); i++)
		{
			file << pos_2d(i,0) << "\t" << pos_2d(i,1) << "\n";
		}
	}
	else
	{
		write_log("error: could not creat output file!");
		exit(EXIT_FAILURE);
	}
	file.close();

	// save coordinates of atoms in 3d space
	file.open(name+".pos_3d.dat", std::ios::trunc);

	if (file.is_open())
	{
		file << std::scientific;
		file << std::showpos;
		file << pos_3d.dim1() << "\t" << pos_3d.dim2() << std::endl;
		for (int i=0; i<pos_3d.dim1(); i++)
		{
			file << pos_3d(i,0) << "\t" << pos_3d(i,1) << "\t" << pos_3d(i,2) << "\n";
		}
	}
	else
	{
		write_log("error: could not creat output file!");
		exit(EXIT_FAILURE);
	}
	file.close();

}

void cnt::electron()
{

	// make the list of 1st nearest neighbor atoms
	nr::mat_int nn_list(2*Nu,3,0); // contains index of the nearest neighbor atom
	nr::mat_int nn_tvec_index(2*Nu,3,0); // contains the index of the cnt unit cell that the nearest neigbor atom is in.
	for (int i=0; i<pos_3d.dim1(); i++)
	{
		int k=0;
		for (int j=0; j<pos_3d.dim1(); j++)
		{
			for (int l=-1; l<=1; l++)
			{
				double dx = pos_3d(i,0) - (pos_3d(j,0) + (double)l*t_vec_3d(0));
				double dy = pos_3d(i,1) - (pos_3d(j,1) + (double)l*t_vec_3d(1));
				double dz = pos_3d(i,2) - (pos_3d(j,2) + (double)l*t_vec_3d(2));
				double dR = sqrt(dx*dx+dy*dy+dz*dz);
				if ( (i!=j) && (dR<(1.4*a_cc)) )
				{
					nn_list(i,k) = j;
					nn_tvec_index(i,k) = l;
					k++;
				}
			}
		}
	}

	nr::mat_complex H(2*Nu, 2*Nu, 0.0);
	nr::mat_complex S(2*Nu, 2*Nu, 0.0);
	nr::vec_doub E(2*Nu,0.0);
	nr::mat_complex C(2*Nu, 2*Nu, 0.0);

	int NK = nk;


	el_energy.assign(2*Nu, NK, 0.0);
	el_psi.assign(2*Nu, 2*Nu, NK, nr::cmplx(0.0));
	
	double t_len = t_vec_3d.norm2();

	// open file to save electron energy bands
	std::ofstream file;
	file.open(name+".electron_energy.dat", std::ios::app);

	for (int n=0; n<NK; n++)
	{
		double wave_vec = double(n-nk/2)*dk_l.norm2();

		H.assign(2*Nu, 2*Nu, 0.0);
		S.assign(2*Nu, 2*Nu, 0.0);

		for (int i=0; i<2*Nu; i++)
		{
			H(i,i) = (e2p,0.e0);
			for (int k=0; k<3; k++)
			{
				int j = nn_list(i,k);
				int l = nn_tvec_index(i,k);

				H(i,j) = H(i,j)+nr::cmplx(t0,0.e0)*exp(nr::cmplx(0.0,wave_vec*double(l)*t_len));
				S(i,j) = S(i,j)+nr::cmplx(s0,0.e0)*exp(nr::cmplx(0.0,wave_vec*double(l)*t_len));
			}
		}

		nr::eig_sym(E, C, H);
		
		// save electron energy bands
		file << std::scientific;
		file << std::showpos;
		file << wave_vec << "\t";
		for (int i=0; i<E.size(); i++)
		{
			file << std::scientific;
			file << std::showpos;
			file << E(i)/t0 << "\t";
		}
		file << "\n";

		// fix the phase of the eigen vectors
		for (int i=0; i<C.dim2(); i++)
		{
			nr::cmplx phi = std::conj(C(0,i))/std::abs(C(0,i));
			for (int j=0; j<C.dim1(); j++)
			{
				C(j,i) = C(j,i)*phi;
			}
		}

		// save electron energy and wavefunction data in member variables
		for (int i=0; i<C.dim2(); i++)
		{
			el_energy(i,n) = E(i);
			for (int j=0; j<C.dim1(); j++)
			{
				el_psi(j,i,n) = C(j,i);
			}
		}

	}
	file.close();

}


void cnt::dielectric()
{


	// calculate v_q
	int Nq = nk;
	nr::mat3d_complex v_q(2*Nu, 2*Nu, Nq, nr::cmplx(0.0));
	nr::vec_doub wave_vec(Nq, 0.0);

	double t_len = t_vec_3d.norm2();
	double factor = nr::sqr(4*constants::pi*constants::eps0*Upp/nr::sqr(constants::q0));
	for(int iq = 0; iq<Nq; iq++)
	{
		wave_vec(iq) = double(iq-Nq/2)*dk_l.norm2();
		for (int i=-int(number_of_cnt_unit_cells/2); i<int(number_of_cnt_unit_cells/2); i++)
		{
			for (int b=0; b<2*Nu; b++)
			{
				for (int bp=0; bp<2*Nu; bp++)
				{
					double dx = pos_3d(b,0) - (pos_3d(bp,0) + double(i)*t_vec_3d(0));
					double dy = pos_3d(b,1) - (pos_3d(bp,1) + double(i)*t_vec_3d(1));
					double dz = pos_3d(b,2) - (pos_3d(bp,2) + double(i)*t_vec_3d(2));
					double dR2 = dx*dx+dy*dy+dz*dz;

					double I = Upp/sqrt(1+factor*dR2);
					v_q(b,bp,iq) += nr::cmplx(1.0/double(number_of_cnt_unit_cells))*exp(nr::cmplx(0.0,wave_vec(iq)*double(i)*t_len))*I;
				}
			}
		}
	}

	// open file to save v_q
	std::ofstream file_vq_real, file_vq_imag;
	file_vq_real.open(name+".vq.real.dat", std::ios::app);
	file_vq_imag.open(name+".vq.imag.dat", std::ios::app);

	file_vq_real << std::scientific;
	file_vq_real << std::showpos;
	file_vq_imag << std::scientific;
	file_vq_imag << std::showpos;
	for (int iq=0; iq<v_q.dim3(); iq++)
	{
		file_vq_real << wave_vec(iq) << "\t";
		file_vq_imag << wave_vec(iq) << "\t";
	}
	file_vq_real << "\n";
	file_vq_imag << "\n";
	for (int b=0; b<v_q.dim1(); b++)
	{
		for (int bp=0; bp<v_q.dim2(); bp++)
		{
			for (int iq=0; iq<v_q.dim3(); iq++)
			{
				file_vq_real << std::real(v_q(b,bp,iq)/constants::eV) << "\t";
				file_vq_imag << std::imag(v_q(b,bp,iq)/constants::eV) << "\t";
			}
			file_vq_real << "\n";
			file_vq_imag << "\n";
		}
	}
	file_vq_real.close();
	file_vq_imag.close();


	// // calculate polarization (PI)
	// nr::vec_doub PI(Nq,0.0);
	// for (int i=0; i<Nq; i++)
	// {
	// 	int iq = i-Nq/2;
	// 	nr::cmplx num1, num2;
	// 	double denom1, denom2;

	// 	for (int ik=0; ik<nk; ik++)
	// 	{
	// 		int ik_p = iq+ik;
	// 		while(ik_p>=nk)	ik_p = ik_p-nk;
	// 		while(ik_p<0) ik_p = ik_p+nk;

	// 		for (int alpha_v=0; alpha_v<Nu; alpha_v++)
	// 		{
	// 			for (int alpha_c=Nu; alpha_c<2*Nu; alpha_c++)
	// 			{
	// 				num1 = nr::cmplx(0.0,0.0);
	// 				num2 = nr::cmplx(0.0,0.0);
	// 				for(int b=0; b<2*Nu; b++)
	// 				{
	// 					num1 += std::conj(el_psi(b,alpha_v,ik)) * el_psi(b,alpha_c,ik_p);
	// 					num2 += std::conj(el_psi(b,alpha_c,ik)) * el_psi(b,alpha_v,ik_p);
	// 				}
	// 				denom1 = el_energy(alpha_c,ik_p)-el_energy(alpha_v,ik);
	// 				denom2 = el_energy(alpha_c,ik)-el_energy(alpha_v,ik_p);

	// 				PI(i) += std::norm(num1)/denom1 + std::norm(num2)/denom2;
	// 			}
	// 		}
	// 	}

	// 	std::cout << "calculating dielectric function: PI(" << i << ") = " << PI(i) << std::endl;
	// }

	// PI = 2.0*PI;

	// // open file to save polarization (PI)
	// std::ofstream file_pi;
	// file_pi.open(name+".pi.dat", std::ios::app);

	// file_pi << std::scientific;
	// file_pi << std::showpos;
	// for (int i=0; i<wave_vec.size(); i++)
	// {
	// 	file_pi << wave_vec(i) << "\t";
	// }
	// file_pi << "\n";
	// for (int i=0; i<PI.size(); i++)
	// {
	// 	file_pi << PI(i) << "\t";
	// }
	// file_pi.close();

	// // calculate dielectric function (epsilon)
	// epsilon.assign(Nq,nr::cmplx(0));
	// for (int i=0; i<epsilon.size(); i++)
	// {
	// 	nr::cmplx sum = nr::cmplx(0);
	// 	for (int b=0; b<2*Nu; b++)
	// 	{
	// 		for (int bp=0; bp<2*Nu; bp++)
	// 		{
	// 			sum += v_q(b,bp,i);
	// 		}
	// 	}
	// 	sum = sum/nr::cmplx((2*Nu)*(2*Nu));
	// 	epsilon(i) = nr::cmplx(1)+sum*nr::cmplx(PI(i));
	// }

	// // open file to save epsilon
	// std::ofstream file_eps;
	// file_eps.open(name+".epsilon.dat", std::ios::app);

	// file_eps << std::scientific;
	// file_eps << std::showpos;
	// for (int i=0; i<wave_vec.size(); i++)
	// {
	// 	file_eps << wave_vec(i) << "\t";
	// }
	// file_eps << "\n";
	// for (int i=0; i<epsilon.size(); i++)
	// {
	// 	file_eps << std::real(epsilon(i)) << "\t";
	// }
	// file_eps << "\n";
	// for (int i=0; i<epsilon.size(); i++)
	// {
	// 	file_eps << std::imag(epsilon(i)) << "\t";
	// }
	// file_eps.close();
}


void cnt::coulomb_int()
{


	// calculate v_q
	int Nq = nk;
	nr::mat3d_complex v_q(2*Nu, 2*Nu, Nq, nr::cmplx(0.0));
	nr::vec_doub wave_vec(Nq, 0.0);

	double t_len = t_vec_3d.norm2();
	double factor = nr::sqr(4*constants::pi*constants::eps0*Upp/nr::sqr(constants::q0));
	for(int iq = 0; iq<Nq; iq++)
	{
		wave_vec(iq) = double(iq)*dk_l.norm2();
		for (int i=-int(number_of_cnt_unit_cells/2); i<int(number_of_cnt_unit_cells/2); i++)
		{
			for (int b=0; b<2*Nu; b++)
			{
				for (int bp=0; bp<2*Nu; bp++)
				{
					double dx = pos_3d(b,0) - (pos_3d(bp,0) + double(i)*t_vec_3d(0));
					double dy = pos_3d(b,1) - (pos_3d(bp,1) + double(i)*t_vec_3d(1));
					double dz = pos_3d(b,2) - (pos_3d(bp,2) + double(i)*t_vec_3d(2));
					double dR2 = dx*dx+dy*dy+dz*dz;

					double I = Upp/sqrt(1+factor*dR2);
					v_q(b,bp,iq) += nr::cmplx(1.0/double(number_of_cnt_unit_cells))*exp(nr::cmplx(0.0,wave_vec(iq)*double(i)*t_len))*I;
				}
			}
		}
	}

	// direct interaction matrix
	nr::mat_complex V_dir(nk,nk,nr::cmplx(0));
	nr::mat_complex V_xch(nk,nk,nr::cmplx(0));
	nr::mat_complex dE(nk,nk,nr::cmplx(0));
	
	// electron bands
	int ac=Nu;
	int acp=Nu;
	int av=Nu-1;
	int avp=Nu-1;

	int iKcm=0; // center of mass wave vector

	for (int ikc=0; ikc<nk; ikc++)
	{
		int ikv = int(ikc+iKcm);
		while(ikv>=nk)	ikv -= nk;
		while(ikv<0) ikv += nk;

		dE(ikc,ikc) = nr::cmplx(el_energy(ac,ikc)-el_energy(av,ikv));

		for (int ikcp=0; ikcp<nk; ikcp++)
		{
			int ikvp = int(ikcp+iKcm);
			while(ikv>=nk)	ikv -= nk;
			while(ikv<0) ikv += nk;

			int iq_dir = ikc-ikcp;
			while(iq_dir>=Nq) iq_dir = iq_dir-Nq;
			while(iq_dir<0) iq_dir = iq_dir+Nq;

			int iq_xch = int(ikc-ikv);
			while(iq_xch>=Nq) iq_xch -= Nq;
			while(iq_xch<0) iq_xch += Nq;

			for (int b=0; b<2*Nu; b++)
			{
				for (int bp=0; bp<2*Nu; bp++)
				{
					V_dir(ikc,ikcp) += std::conj(el_psi(b,ac,ikc))*el_psi(b,acp,ikcp)*el_psi(bp,av,ikv)*std::conj(el_psi(bp,avp,ikvp))*v_q(b,bp,iq_dir);
					V_xch(ikc,ikcp) += std::conj(el_psi(b,ac,ikc))*el_psi(b,av,ikv)*el_psi(bp,acp,ikcp)*std::conj(el_psi(bp,avp,ikvp))*v_q(b,bp,iq_xch);
				}
			}
		}
	}


	// open file to save V_dir
	std::ofstream file_V_dir_real, file_V_dir_imag;
	file_V_dir_real.open(name+".v_dir.real.dat", std::ios::app);
	file_V_dir_imag.open(name+".v_dir.imag.dat", std::ios::app);

	file_V_dir_real << std::scientific << std::showpos;
	file_V_dir_imag << std::scientific << std::showpos;
	
	file_V_dir_real << double(0.0) << "\t";
	file_V_dir_imag << double(0.0) << "\t";
	for (int ikcp=0; ikcp<V_dir.dim2(); ikcp++)
	{
		double wave_vec = double(ikcp)*dk_l.norm2();
		file_V_dir_real << wave_vec << "\t";
		file_V_dir_imag << wave_vec << "\t";
	}
	file_V_dir_real << "\n";
	file_V_dir_imag << "\n";

	for (int ikc=0; ikc<V_dir.dim1(); ikc++)
	{
		double wave_vec = double(ikc)*dk_l.norm2();
		file_V_dir_real << wave_vec << "\t";
		file_V_dir_imag << wave_vec << "\t";

		for (int ikcp=0; ikcp<V_dir.dim2(); ikcp++)
		{
			file_V_dir_real << std::real(V_dir(ikc,ikcp)/constants::eV) << "\t";
			file_V_dir_imag << std::imag(V_dir(ikc,ikcp)/constants::eV) << "\t";
		}
		file_V_dir_real << "\n";
		file_V_dir_imag << "\n";
	}
	file_V_dir_real.close();
	file_V_dir_imag.close();

	// open file to save V_xch
	std::ofstream file_V_xch_real, file_V_xch_imag;
	file_V_xch_real.open(name+".v_xch.real.dat", std::ios::app);
	file_V_xch_imag.open(name+".v_xch.imag.dat", std::ios::app);

	file_V_xch_real << std::scientific << std::showpos;
	file_V_xch_imag << std::scientific << std::showpos;
	
	file_V_xch_real << double(0.0) << "\t";
	file_V_xch_imag << double(0.0) << "\t";
	for (int ikcp=0; ikcp<V_xch.dim2(); ikcp++)
	{
		double wave_vec = double(ikcp)*dk_l.norm2();
		file_V_xch_real << wave_vec << "\t";
		file_V_xch_imag << wave_vec << "\t";
	}
	file_V_xch_real << "\n";
	file_V_xch_imag << "\n";

	for (int ikc=0; ikc<V_xch.dim1(); ikc++)
	{
		double wave_vec = double(ikc)*dk_l.norm2();
		file_V_xch_real << wave_vec << "\t";
		file_V_xch_imag << wave_vec << "\t";

		for (int ikcp=0; ikcp<V_xch.dim2(); ikcp++)
		{
			file_V_xch_real << std::real(V_xch(ikc,ikcp)/constants::eV) << "\t";
			file_V_xch_imag << std::imag(V_xch(ikc,ikcp)/constants::eV) << "\t";
		}
		file_V_xch_real << "\n";
		file_V_xch_imag << "\n";
	}
	file_V_xch_real.close();
	file_V_xch_imag.close();


	// solve for exciton energy
	nr::mat_complex kernel(dE-V_dir);
	// nr::mat_complex kernel(dE+nr::cmplx(2)*V_xch-V_dir);
	
	// nr::mat_complex ex_psi(nk,nk);
	// nr::vec_doub ex_energy(nk);
	// nr::eig_sym(ex_energy, ex_psi, kernel);

	nr::mat_complex ex_psi(nk,nk,nr::cmplx(0));
	nr::vec_doub ex_energy(nk,0.0);
	nr::eig_sym_selected(ex_energy, ex_psi, kernel, 20);


	// open file to save exciton energies
	std::ofstream file_ex_energy;
	file_ex_energy.open(name+".ex_energy.dat", std::ios::app);
	file_ex_energy << std::scientific << std::showpos;
	for (int i=0; i<ex_energy.size(); i++)
	{
		file_ex_energy << ex_energy(i)/constants::eV << "\t" << std::real(dE(i,i))/constants::eV << "\n";
	}
	file_ex_energy.close();


	// open file to save exciton energies
	std::ofstream file_ex_psi;
	file_ex_psi.open(name+".ex_psi.dat", std::ios::app);
	file_ex_psi << std::scientific << std::showpos;
	for (int i=0; i<ex_psi.dim1(); i++)
	{
		for (int j=0; j<ex_psi.dim2(); j++)
		{
			file_ex_psi << std::real(ex_psi(i,j)) << "\t";
		}
		file_ex_psi << "\n";
	}
	file_ex_psi.close();

}