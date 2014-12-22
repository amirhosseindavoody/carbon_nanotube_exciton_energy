subroutine fnSelfEnergy()
  use comparams
  implicit none
  
  integer :: ik, mu_k, iq, mu_q, ikq, mu_kq
  integer :: tmpi,i,u,s
  real*8, dimension(2) :: k,kq
  real*8, dimension(2) :: E_kq
  complex*16, dimension(2) :: Cc_kq,Cv_kq
  
  ! initialize variables.**********************************************************************************************
  allocate(Ek(2,ik_low:ik_high,2))
  allocate(Sk(2,ik_low:ik_high,2))
  allocate(Cc(2,ik_low:ik_high,2))
  allocate(Cv(2,ik_low:ik_high,2))
  
  do ik=ik_low,ik_high
    Sk(1,ik,1)=0.d0
    Sk(1,ik,2)=0.d0
    Sk(2,ik,1)=0.d0
    Sk(2,ik,2)=0.d0
  enddo
  
  ! calculate the self energies and the tight-binding energies and coefficients.***************************************
  do ik=ik_low,ik_high
    print *,'self-energy: ik=',ik
    mu_k=min_sub(i_sub)
    k=dble(mu_k)*K1+dble(ik)*dk*K2
    call fnGrapheneEnergy(Ek(1,ik,:),Cc(1,ik,:),Cv(1,ik,:),k)
    
    !calculate self energy of the first subband
    do mu_q=1-Nu/2,Nu/2
      do iq=ikc_min,ikc_max
        kq=dble(mu_k+mu_q)*K1+dble(ik+iq)*dk*K2
        call fnGrapheneEnergy(E_kq,Cc_kq,Cv_kq,kq) 
        Sk(1,ik,1)=Sk(1,ik,1)-dreal(conjg(Cc(1,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,1)*Cc(1,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cc(1,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,2)*Cc(1,ik,2)*conjg(Cv_kq(2))+ &
                                    conjg(Cc(1,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,1)*Cc(1,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cc(1,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,2)*Cc(1,ik,2)*conjg(Cv_kq(2)))/(kappa*eps_q(-mu_q,-iq))
        
        Sk(1,ik,2)=Sk(1,ik,2)-dreal(conjg(Cv(1,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,1)*Cv(1,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cv(1,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,2)*Cv(1,ik,2)*conjg(Cv_kq(2))+ &
                                    conjg(Cv(1,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,1)*Cv(1,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cv(1,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,2)*Cv(1,ik,2)*conjg(Cv_kq(2)))/(kappa*eps_q(-mu_q,-iq))
      enddo
    enddo
    
    !calculate the self energy of the second subband
    mu_k=-min_sub(i_sub)
    k=dble(mu_k)*K1+dble(ik)*dk*K2
    call fnGrapheneEnergy(Ek(2,ik,:),Cc(2,ik,:),Cv(2,ik,:),k)
    
    do mu_q=1-Nu/2,Nu/2
      do iq=ikc_min,ikc_max
        kq=dble(mu_k+mu_q)*K1+dble(ik+iq)*dk*K2
        call fnGrapheneEnergy(E_kq,Cc_kq,Cv_kq,kq) 
        Sk(2,ik,1)=Sk(2,ik,1)-dreal(conjg(Cc(2,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,1)*Cc(2,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cc(2,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,2)*Cc(2,ik,2)*conjg(Cv_kq(2))+ &
                                    conjg(Cc(2,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,1)*Cc(2,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cc(2,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,2)*Cc(2,ik,2)*conjg(Cv_kq(2)))/(kappa*eps_q(-mu_q,-iq))
        
        Sk(2,ik,2)=Sk(2,ik,2)-dreal(conjg(Cv(2,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,1)*Cv(2,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cv(2,ik,1))*Cv_kq(1)*v_FT(-mu_q,-iq,1,2)*Cv(2,ik,2)*conjg(Cv_kq(2))+ &
                                    conjg(Cv(2,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,1)*Cv(2,ik,1)*conjg(Cv_kq(1))+ &
                                    conjg(Cv(2,ik,2))*Cv_kq(2)*v_FT(-mu_q,-iq,2,2)*Cv(2,ik,2)*conjg(Cv_kq(2)))/(kappa*eps_q(-mu_q,-iq))
      enddo
    enddo
  enddo
  
  !save the results.***************************************************************************************************
  do ik=ik_low,ik_high
    write(fh11,10, advance='no') Sk(1,ik,1) 
    write(fh12,10, advance='no') Sk(1,ik,2) 
  enddo
  
  write(fh11,10)
  write(fh12,10)
  
  do ik=ik_low,ik_high
    write(fh11,10, advance='no') Sk(2,ik,1) 
    write(fh12,10, advance='no') Sk(2,ik,2) 
  enddo
  
10 FORMAT (E16.8)
  
  
  return
end