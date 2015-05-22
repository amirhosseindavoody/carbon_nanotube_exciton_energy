!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of physical variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physical_constant_mod
  implicit none
  private
  
  !Mathematical constants
  real*8, public, parameter :: pi=3.141592d0
  complex*16, public, parameter :: i1=(0.d0,1.d0)
  
  !Physical constants
  real*8, public, parameter :: eV=1.6d-19 ![Joules]
  real*8, public, parameter :: hb=6.5d-16*eV ![eV.s]
  real*8, public, parameter :: kb=1.3865d-23 ![J/K]
  
  real*8, public, parameter :: a_cc=1.42d-10 !carbon-carbon distance [meters] 
  real*8, public, parameter :: a_l=dsqrt(3.d0)*a_cc !lattice constants
  
  real*8, public, parameter :: e2p = 0.d0 !tight binding constants
  real*8, public, parameter :: t0 = 2.7d0*eV  !tight binding constants
  real*8, public, parameter :: s0 = 0.d0 !tight binding constants
  
  real*8, public, parameter :: Upp = 11.3d0*eV !constant used in the Ohno potential
  real*8, public, parameter :: eps0 = 8.85d-12 !permittivity of free space
  real*8, public, parameter :: q0 = 1.6d-19 !charge of electron
  
end module physical_constant_mod