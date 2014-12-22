!*******************************************************************************
! This subroutines sets the physical constants of the simulation
!*******************************************************************************
subroutine fnPhysConst()
  use comparams
  implicit none
  
  eV=1.6d-19 ![Joules]
  hb=6.5d-16*eV ![eV.s]
  a_cc=1.42d-10 ![meters]
  a_l=dsqrt(3.d0)*a_cc
  
  e2p=0 ![Joules]
  t0=2.7d0*eV
  s0=0.d0
  
  Upp=11.3d0*eV
  eps0=8.85d-12
  q0=1.6d-19
  
  kappa=2.d0
  
  E_th=E_th*eV
  
  
  write(fh1,10) "TIGHTBINDING PARAMETERS"
  write(fh1,11) "e2p [eV]=",e2p/eV
  write(fh1,11) "t0 [eV]=",t0/eV
  write(fh1,11) "s0 [eV]=",s0
  
10 FORMAT (A100)
11 FORMAT (A10,E16.8)   
12 FORMAT (A10,I3) 
  return
end