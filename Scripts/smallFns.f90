!**********************************************************************************************************************
! This subroutine calculates the greatest common divisor of the arguments na and nb
!**********************************************************************************************************************
subroutine gcd(ngcd,na,nb)
  integer :: na,nb,ngcd
  integer :: ia,ib,itemp
    
  ia=na
  ib=nb
do while (ib .ne. 0)
      itemp=ia
      ia=ib
      ib=mod(itemp,ib)
  end do
  ngcd=ia
  return
end

!**********************************************************************************************************************
! This subroutine saves the some of the simulation parameters
!**********************************************************************************************************************
subroutine fnSaveSimInfo()
  use comparams
  implicit none

  write(fh1,10) "SIMULATION PARAMETERS"
  write(fh1,12) "n_ch=",n_ch
  write(fh1,12) "m_ch=",m_ch
  write(fh1,12) "nkg=",nkg
  write(fh1,12) "nr=",nr
  write(fh1,14) "E_th=",E_th
  write(fh1,14) "Kcm_max=",Kcm_max
  write(fh1,13) "flg_dielectric=",flg_dielectric
  write(fh1,12) "i_sub=",i_sub
  
10 FORMAT (A100)
11 FORMAT (A10,E16.8,E16.8)   
12 FORMAT (A10,I5) 
13 FORMAT (A15,L1) 
14 FORMAT (A10,E16.8)  

  return
end
    
!**********************************************************************************************************************
! This subroutine saves the limits of different types of k-vectors in the simulated exciton energy code.
!**********************************************************************************************************************
subroutine fnSaveMisc()
  use comparams
  implicit none

  write(fh19,10) min_sub(i_sub)
  write(fh19,10) ik_max
  write(fh19,10) iKcm_max
  write(fh19,10) ikr_high
  write(fh19,10) ik_high
  write(fh19,10) iq_max
  write(fh19,10) nX
  write(fh19,11) dk
  
10 FORMAT (I4.4) 
11 FORMAT (E16.8)

  return
end