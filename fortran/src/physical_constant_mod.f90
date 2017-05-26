!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of physical variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physical_constant_mod
  implicit none
  private
  
  !Mathematical constants
  double precision, public, parameter :: pi=3.141592d0
  double complex, public, parameter :: i1=(0.d0,1.d0)
  
  !Physical constants
  double precision, public, parameter :: eV=1.6d-19 ![Joules]
  double precision, public, parameter :: hb=6.5d-16*eV ![eV.s]
  double precision, public, parameter :: kb=1.3865d-23 ![J/K]
  
  double precision, public, parameter :: a_cc=1.42d-10 !carbon-carbon distance [meters] 
  double precision, public, parameter :: a_l=dsqrt(3.d0)*a_cc !lattice constants
  
  double precision, public, parameter :: e2p = 0.d0 !tight binding constants
  double precision, public, parameter :: t0 = 2.7d0*eV  !tight binding constants
  double precision, public, parameter :: s0 = 0.d0 !tight binding constants
  
  double precision, public, parameter :: Upp = 11.3d0*eV !constant used in the Ohno potential
  double precision, public, parameter :: eps0 = 8.85d-12 !permittivity of free space
  double precision, public, parameter :: q0 = 1.6d-19 !charge of electron

  double precision, public, parameter :: kappa = 2.0d0 !dielectric constant due to surrounding media
  
end module physical_constant_mod