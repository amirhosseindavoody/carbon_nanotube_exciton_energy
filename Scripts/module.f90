!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of global variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module comparams
  use dfport ! Standard library for rand
  INCLUDE 'link_fnl_static.h'
  implicit none
  
  real*8, parameter :: pi=3.141592d0
  complex*16, parameter :: i1=(0.d0,1.d0)
  
  !Input parameters
  integer :: n_ch,m_ch !chiral vector parameters
  integer :: nkg,nr,i_sub !reciprocal space mesh size in graphene and number of CNT unit cells in real space
  real*8 :: E_th !threshold energy
  real*8 :: Kcm_max
  
  !Flags
  logical :: flg_dielectric
  
  !Physical constants
  real*8 :: eV,hb
  real*8 :: a_cc,a_l !lattice constants
  real*8 :: e2p,t0,s0 !tight binding constants
  real*8 :: Upp !constant used in the Ohno potential
  real*8 :: eps0 !permittivity of free space
  real*8 :: q0 !charge of electron
  real*8 :: kappa !dielectric constant due to core electrons in CNT
  
  !Geometrical properties
  real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
  real*8 :: len_ch,radius
  integer :: dR,Nu,MC
  integer :: t1,t2
  real*8, dimension(:,:), allocatable :: posA,posB,posAA,posBB,posAB,posBA
  
  !Reciprocal lattice properties
  real*8 :: dk
  real*8, dimension(2) :: K1, K2
  
  !CNT band structure properties
  integer, dimension(:), allocatable :: min_sub
  integer :: ikc_max, ikc_min, ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
  
  !Dielectric function
  real*8, dimension(:,:), allocatable :: eps_q
  complex*16, dimension(:,:,:,:), allocatable :: v_FT
  
  !CNT self energy and tight binding coefficients
  real*8, dimension(:,:,:), allocatable:: Ek,Sk
  complex*16, dimension(:,:,:), allocatable:: Cc,Cv
  
  integer :: nX ! this is the number of exciton bands below free quasi-particle bands
  
  
  !These are the unit numbers for the input/output files in the program
  integer :: fh1,fh2,fh3,fh4,fh5,fh6,fh7,fh8,fh9,fh10,fh11,fh12,fh13,fh14,fh15,fh16,fh17,fh18,fh19
  
    end module comparams
    
    
! Explanation of variables
    ! Ek(mu,k,n) stores the tight-binding energy of the band "mu" with wavevector "k" in conduction band (n=1)
    !               or the valence band (n=2).
    ! v_FT(mu,q,n,m) stores the Fourier transform of the Coulomb potential at the wavevector determined
    !               by band index "mu" and wavenumber "q" for atoms of type A (n or m = 1) or type B (n or m = 2)
    

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of exciton energy and wavefunction variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module exciton_prop
  implicit none
  
  real*8, dimension(:), allocatable :: Ex_A1, Ex0_A2, Ex1_A2
  complex*16, dimension(:,:), allocatable :: Psi_A1, Psi0_A2, Psi1_A2
  
end module exciton_prop    