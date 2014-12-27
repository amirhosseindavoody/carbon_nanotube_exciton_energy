subroutine fnExcitonDispersion()
  use comparams
  use exciton_prop
  implicit none
  
  integer :: iKcm, ikr, ikpr
  real*8 :: Ef_min
  
  allocate(Psi_A1(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Psi0_A2(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Psi1_A2(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Ex_A1(ikr_low:ikr_high))
  allocate(Ex0_A2(ikr_low:ikr_high))
  allocate(Ex1_A2(ikr_low:ikr_high))
  
  nX=0
  iKcm=0
  call fnExcitonEnergy(Ef_min, iKcm)
  
  do ikr=ikr_low,ikr_high
    if (Ef_min .ge. Ex0_A2(ikr)) then
      nX=nX+1
    endif
  enddo
  
  
  
  do iKcm=iKcm_min,iKcm_max
    print *,"Exciton Disp. --> iKcm=", iKcm
    
    call fnExcitonEnergy(Ef_min, iKcm)
    
    ! save exciton energy and wavefunction
    do ikr=ikr_low,(ikr_low+nX-1)
      write(fh13,10, advance='no') Ex_A1(ikr)
      write(fh14,10, advance='no') Ex0_A2(ikr)
      write(fh15,10, advance='no') Ex1_A2(ikr)
      do ikpr=ikr_low,ikr_high
        write(fh16,11, advance='no') Psi_A1(ikpr,ikr)
        write(fh17,11, advance='no') Psi0_A2(ikpr,ikr)
        write(fh18,11, advance='no') Psi1_A2(ikpr,ikr)
      enddo
    enddo
    write(fh13,10)
    write(fh14,10)
    write(fh15,10)
    write(fh16,10)
    write(fh17,10)
    write(fh18,10)
  
  enddo
  
  deallocate(Psi_A1)
  deallocate(Psi0_A2)
  deallocate(Psi1_A2)
  deallocate(Ex_A1)
  deallocate(Ex0_A2)
  deallocate(Ex1_A2)

10 FORMAT (E16.8)  
11 FORMAT (E16.8,E16.8) 
  
  return
end