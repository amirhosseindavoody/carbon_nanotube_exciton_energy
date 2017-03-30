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
// #include "write_log.h"
#include "nr3.h"

cnt::cnt(const std::string &in_name, const int in_n, const int in_m, const int in_length)
{
	name = in_name;
	n = in_n;
	m = in_m;
	length_in_cnt_unit_cell = in_length;

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

	std::cout << "a1 = " << a1[0] << " , " << a1[1] << std::endl;
	std::cout << "a2 = " << a2[0] << " , " << a2[1] << std::endl;
	std::cout << "aCC_vec = " << aCC_vec[0] << " , " << aCC_vec[1] << std::endl;

	
	// calculate chirality and translational vectors of CNT unit cell
	ch_vec = (double)n * a1 + (double)m * a2;
	ch_len = a_l*sqrt(pow((double)n,2)+pow((double)m,2)+(double)n*m);
	radius = ch_len/2.0/constants::pi;

	int dR = nr::gcd(2*n+m,n+2*m);
	int t1 = +(2*m+n)/dR;
	int t2 = -(2*n+m)/dR;
	t_vec = (double)t1*a1 + (double)t2*a2;

	Nu = 2*(pow(n,2)+pow(m,2)+n*m)/dR;

	// // rotate basis vectors so that ch_vec is along the x_axis
	// double cos_theta = ch_vec[0]/ublas::norm_2(ch_vec);
	// double sin_theta = ch_vec[1]/ublas::norm_2(ch_vec);
	// ublas::matrix <double> rot(2,2); // rotation matrix
	// rot(0,0) = cos_theta;
	// rot(1,0) = -sin_theta;
	// rot(0,1) = sin_theta;
	// rot(1,1) = cos_theta;

	// ch_vec = ublas::prod(rot,ch_vec);
	// t_vec = ublas::prod(rot,t_vec);
	// a1 = ublas::prod(rot,a1);
	// a2 = ublas::prod(rot,a2);
	// b1 = ublas::prod(rot,b1);
	// b2 = ublas::prod(rot,b2);
	// aCC_vec = ublas::prod(rot,aCC_vec);

	// // calculate reciprocal lattice of CNT
	// K1 = (-t2*b1 + t1*b2)/((double)Nu);
	// K2 = ((double)m*b1-(double)n*b2)/((double)Nu);
	// dk_l = K1/((double)length_in_cnt_unit_cell);

	// // calculate positions of atoms in the cnt unit cell
	// pos_a.resize(Nu,2);
	// pos_b.resize(Nu,2);

	// int k =0;

	// for (int i=0; i<=t1+n; i++)
	// {
	// 	for (int j=t2; j<=m; j++)
	// 	{
	// 		bool flag1 = (double)(t2*i)/(double)t1 <= (double)j;
	// 		bool flag2 = (double)(m*i)/(double)n >= (double)j;
	// 		bool flag3 = (double)(t2*(i-n))/(double)(t1) > (double)(j-m);
	// 		bool flag4 = (double)(m*(i-t1))/(double)(n) < (double)(j-t2);

	// 		if(flag1 && flag2 && flag3 && flag4)
	// 		{
	// 			pos_a(k,0) = (double)i*a1(0) + (double)j*a2(0);
	// 			pos_a(k,1) = (double)i*a1(1) + (double)j*a2(1);
	// 			pos_b(k,0) = pos_a(k,0)+aCC_vec(0);
	// 			pos_b(k,1) = pos_a(k,1)+aCC_vec(1);

	// 			if(pos_a(k,0) > ch_vec(0))	pos_a(k,0) -= ch_vec(0);
	// 			if(pos_a(k,0) < 0.0)	pos_a(k,0) += ch_vec(0);
	// 			if(pos_a(k,1) > ch_vec(1))	pos_a(k,1) -= ch_vec(1);
	// 			if(pos_a(k,1) < 0.0)	pos_a(k,1) += ch_vec(1);

	// 			if(pos_b(k,0) > ch_vec(0))	pos_b(k,0) -= ch_vec(0);
	// 			if(pos_b(k,0) < 0.0)	pos_b(k,0) += ch_vec(0);
	// 			if(pos_b(k,1) > ch_vec(1))	pos_b(k,1) -= ch_vec(1);
	// 			if(pos_b(k,1) < 0.0)	pos_b(k,1) += ch_vec(1);
				
	// 			k++;
	// 		}
	// 	}
	// }

	// if (k != Nu) 
	// {
	// 	cout << "error in finding position of atoms in cnt unit cell!!!" << endl;
	// 	cout << "Nu = " << Nu << "  ,  k = " << k << endl;
	// 	exit(1);
	// }

	// // // calculate distances between atoms in a warped cnt unit cell.
	// // pos_aa.resize(Nu,2);
	// // pos_ab.resize(Nu,2);
	// // pos_ba.resize(Nu,2);
	// // pos_bb.resize(Nu,2);

	// // for (int i=0; i<Nu; i++)
	// // {
	// // 	pos_aa(i,0) = pos_a(i,0)-pos_a(0,0);
	// // 	pos_aa(i,1) = pos_a(i,1)-pos_a(0,1);

	// // 	pos_ab(i,0) = pos_a(i,0)-pos_b(0,0);
	// // 	pos_ab(i,1) = pos_a(i,1)-pos_b(0,1);

	// // 	pos_ba(i,0) = pos_b(i,0)-pos_a(0,0);
	// // 	pos_ba(i,1) = pos_b(i,1)-pos_a(0,1);

	// // 	pos_bb(i,0) = pos_b(i,0)-pos_b(0,0);
	// // 	pos_bb(i,1) = pos_b(i,1)-pos_b(0,1);

	// // 	if(pos_aa(i,0) > ch_vec(0)/2)	pos_aa(i,0) -= ch_vec(0);
	// // 	if(pos_ab(i,0) > ch_vec(0)/2)	pos_ab(i,0) -= ch_vec(0);
	// // 	if(pos_ba(i,0) > ch_vec(0)/2)	pos_ba(i,0) -= ch_vec(0);
	// // 	if(pos_bb(i,0) > ch_vec(0)/2)	pos_bb(i,0) -= ch_vec(0);
	// // }

	// // put position of all atoms in a single variable in 2d space(unrolled graphene sheet)
	// pos_2d.resize(2*Nu,2);
	// for (int i=0; i<Nu; i++)
	// {
	// 	pos_2d(i,0) = pos_a(i,0);
	// 	pos_2d(i+Nu,0) = pos_b(i,0);
	// 	pos_2d(i,1) = pos_a(i,1);
	// 	pos_2d(i+Nu,1) = pos_b(i,1);
	// }

	// // calculate position of all atoms in the 3d space (rolled graphene sheet)
	// pos_3d.resize(2*Nu,3);
	// for (int i=0; i<pos_3d.size1(); i++)
	// {
	// 	pos_3d(i,0) = radius*cos(pos_2d(i,0)/radius);
	// 	pos_3d(i,1) = pos_2d(i,1);
	// 	pos_3d(i,2) = radius*sin(pos_2d(i,0)/radius);
	// }

	// //make 3d t_vec
	// t_vec_3d.resize(3);
	// t_vec_3d(0) = 0;
	// t_vec_3d(1) = t_vec(1);
	// t_vec_3d(2) = 0;


	// // save coordinates of atoms in 2d space
	// ofstream file;
	// file.open(name+".pos_2d.dat", ios::trunc);

	// if (file.is_open())
	// {
	// 	file << std::scientific;
	// 	file << std::showpos;
	// 	file << pos_2d.size1() << "\t" << pos_2d.size2() << endl;
	// 	for (int i=0; i<pos_2d.size1(); i++)
	// 	{
	// 		file << pos_2d(i,0) << "\t" << pos_2d(i,1) << "\n";
	// 	}
	// }
	// else
	// {
	// 	write_log("error: could not creat output file!");
	// 	exit(EXIT_FAILURE);
	// }
	// file.close();

	// // save coordinates of atoms in 3d space
	// file.open(name+".pos_3d.dat", ios::trunc);

	// if (file.is_open())
	// {
	// 	file << std::scientific;
	// 	file << std::showpos;
	// 	file << pos_3d.size1() << "\t" << pos_3d.size2() << endl;
	// 	for (int i=0; i<pos_3d.size1(); i++)
	// 	{
	// 		file << pos_3d(i,0) << "\t" << pos_3d(i,1) << "\t" << pos_3d(i,2) << "\n";
	// 	}
	// }
	// else
	// {
	// 	write_log("error: could not creat output file!");
	// 	exit(EXIT_FAILURE);
	// }
	// file.close();

}

void cnt::electron()
{
	// namespace ublas = boost::numeric::ublas;

	// // make the list of 1st nearest neighbor atoms
	// ublas::matrix<int> nn_list(2*Nu,3); // contains index of the nearest neighbor atom
	// ublas::matrix<int> nn_tvec_index(2*Nu,3); // contains the index of the cnt unit cell that the nearest neigbor atom is in.
	// for (int i=0; i<pos_3d.size1(); i++)
	// {
	// 	int k=0;
	// 	for (int j=0; j<pos_3d.size1(); j++)
	// 	{
	// 		for (int l=-1; l<=1; l++)
	// 		{
	// 			double dx = pos_3d(i,0) - pos_3d(j,0) + (double)l*t_vec_3d(0);
	// 			double dy = pos_3d(i,1) - pos_3d(j,1) + (double)l*t_vec_3d(1);
	// 			double dz = pos_3d(i,2) - pos_3d(j,2) + (double)l*t_vec_3d(2);
	// 			double dR = sqrt(dx*dx+dy*dy+dz*dz);
	// 			if ( (i!=j) && (dR<(1.4*a_cc)) )
	// 			{
	// 				nn_list(i,k) = j;
	// 				nn_tvec_index(i,k) = l;
	// 				k++;
	// 			}
	// 		}
	// 	}
	// }

	// // ublas::vector<std::complex<double>> C(2*Nu);
	// // ublas::matrix<std::complex<double>> H(2*Nu,2*Nu);
	// // ublas::matrix<std::complex<double>> S(2*Nu,2*Nu);
	
	// // double t_len = ublas::norm_2(t_vec_3d);

	// // C.clear();
	// // H.clear();
	// // S.clear();


	// // for (int i=0; i<2*Nu; i++)
	// // {
	// // 	H(i,i) = (e2p,0.e0);
	// // 	for (int k=0; k<3; k++)
	// // 	{
	// // 		int j = nn_list(i,k);
	// // 		int l = nn_tvec_index(i,k);

	// // 		H(i,j) = (t0,0.e0)*exp(1i*((double)l*t_len,0.e0));
	// // 		S(i,j) = (s0,0.e0)*exp(1i*((double)l*t_len,0.e0));
	// // 	}
	// // }

	// // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> m(2,2);
	// // m.resize(2,2);

	// // my_math::c_vector c(2*Nu);
	// // my_math::c_matrix H(2*Nu, 2*Nu);
	// // my_math::c_matrix S(2*Nu, 2*Nu);
	
	// // // double t_len = ublas::norm_2(t_vec_3d);

	// // // C.clear();
	// // // H.clear();
	// // // S.clear();

	// // std::cout << "C rows = " << c.size1() << std::endl;
	// // std::cout << "C size = " << H.size1() << "x" << H.size2() << std::endl;


	// // for (int i=0; i<2*Nu; i++)
	// // {
	// // 	H(i,i) = (e2p,0.e0);
	// // 	for (int k=0; k<3; k++)
	// // 	{
	// // 		int j = nn_list(i,k);
	// // 		int l = nn_tvec_index(i,k);

	// // 		H(i,j) = (t0,0.e0)*exp(1i*((double)l*t_len,0.e0));
	// // 		S(i,j) = (s0,0.e0)*exp(1i*((double)l*t_len,0.e0));
	// // 	}
	// // }


}