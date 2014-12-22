subroutine fnExcitonEnergy(Ef_min, iKcm)
  use comparams, only: Cc, Cv, Ek, Sk, v_FT, eps_q, kappa,min_sub,i_sub,ikr_high,ikr_low,dk,K1,K2
  use exciton_prop
  use IMSL_libraries
  implicit none
  
  integer :: nkr
  integer :: iKcm,mu_cm,mu_kr
  integer :: ikr, ikc, ikv, ikpr, ikpc, ikpv,i
  real :: error
  real*8 :: tmpr, Ef_min
  real*8, dimension(2) :: Kcm
  real*8, dimension(:), allocatable :: Ef
  complex*16, dimension(:), allocatable :: Psi_tmp
  complex*16, dimension(:,:), allocatable :: Kd11, Kd12, Kx11, Ke11, Kernel, Ktmp

  mu_cm=0
  mu_kr=min_sub(i_sub)
  nkr=ikr_high-ikr_low+1
  

  
  allocate(Kd11(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Kd12(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Kx11(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Ke11(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Kernel(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Ktmp(ikr_low:ikr_high,ikr_low:ikr_high))
  allocate(Ef(ikr_low:ikr_high))
  allocate(Psi_tmp(ikr_low:ikr_high))
  
  Kcm=dble(mu_cm)*K1+dble(iKcm)*dk*K2
    
  Ke11=0.d0*Ke11
  Kd11=0.d0*Kd11
  Kd12=0.d0*Kd12
  Kx11=0.d0*Kx11
  
  ! calculate the kernel matrices
  do ikr=ikr_low,ikr_high
    ikc=ikr+iKcm
    ikv=ikr-iKcm
    Ke11(ikr,ikr)=Ek(1,ikc,1)-Ek(1,ikv,2)+Sk(1,ikc,1)-Sk(1,ikv,2)
    Ef(ikr)=Ek(1,ikc,1)-Ek(1,ikv,2)+Sk(1,ikc,1)-Sk(1,ikv,2)
    do ikpr=ikr,ikr_high
      ikpc=ikpr+iKcm
      ikpv=ikpr-iKcm
      Kd11(ikr,ikpr)=(conjg(Cc(1,ikc,1))*Cc(1,ikpc,1)*v_FT(0,ikr-ikpr,1,1)*Cv(1,ikv,1)*conjg(Cv(1,ikpv,1))+ &
                      conjg(Cc(1,ikc,1))*Cc(1,ikpc,1)*v_FT(0,ikr-ikpr,1,2)*Cv(1,ikv,2)*conjg(Cv(1,ikpv,2))+ &
                      conjg(Cc(1,ikc,2))*Cc(1,ikpc,2)*v_FT(0,ikr-ikpr,2,1)*Cv(1,ikv,1)*conjg(Cv(1,ikpv,1))+ &
                      conjg(Cc(1,ikc,2))*Cc(1,ikpc,2)*v_FT(0,ikr-ikpr,2,2)*Cv(1,ikv,2)*conjg(Cv(1,ikpv,2)))/(kappa*eps_q(0,ikr-ikpr))
        
      Kx11(ikr,ikpr)=(conjg(Cc(1,ikc,1))*Cv(1,ikv,1)*v_FT(0,2*iKcm,1,1)*Cc(1,ikpc,1)*conjg(Cv(1,ikpv,1))+ &
                      conjg(Cc(1,ikc,1))*Cv(1,ikv,1)*v_FT(0,2*iKcm,1,2)*Cc(1,ikpc,2)*conjg(Cv(1,ikpv,2))+ &
                      conjg(Cc(1,ikc,2))*Cv(1,ikv,2)*v_FT(0,2*iKcm,2,1)*Cc(1,ikpc,1)*conjg(Cv(1,ikpv,1))+ &
                      conjg(Cc(1,ikc,2))*Cv(1,ikv,2)*v_FT(0,2*iKcm,2,2)*Cc(1,ikpc,2)*conjg(Cv(1,ikpv,2)))
    enddo
    
    do ikpr=ikr_high,-ikr,-1
      ikpc=ikpr+iKcm
      ikpv=ikpr-iKcm
      Kd12(ikr,-ikpr)=(conjg(Cc(1,ikc,1))*Cc(2,ikpc,1)*v_FT(2*mu_kr,ikr-ikpr,1,1)*Cv(1,ikv,1)*conjg(Cv(2,ikpv,1))+ &
                       conjg(Cc(1,ikc,1))*Cc(2,ikpc,1)*v_FT(2*mu_kr,ikr-ikpr,1,2)*Cv(1,ikv,2)*conjg(Cv(2,ikpv,2))+ &
                       conjg(Cc(1,ikc,2))*Cc(2,ikpc,2)*v_FT(2*mu_kr,ikr-ikpr,2,1)*Cv(1,ikv,1)*conjg(Cv(2,ikpv,1))+ &
                       conjg(Cc(1,ikc,2))*Cc(2,ikpc,2)*v_FT(2*mu_kr,ikr-ikpr,2,2)*Cv(1,ikv,2)*conjg(Cv(2,ikpv,2)))/(kappa*eps_q(2*mu_kr,ikr-ikpr))
    enddo
  enddo
  
  ! when running the code in release mode there the next few lines generated a stack overflow error which forced me to use a dummy variable to resolve this issue.
  Ktmp=conjg(transpose(Kd11))
  Kd11=Kd11+Ktmp
  Ktmp=conjg(transpose(Kx11))
  Kx11=Kx11+Ktmp
  Ktmp=conjg(transpose(Kd12))
  Kd12=Kd12+Ktmp
  
  do ikr=ikr_low,ikr_high
    Kd11(ikr,ikr)=Kd11(ikr,ikr)/2.d0
    Kx11(ikr,ikr)=Kx11(ikr,ikr)/2.d0
    Kd12(ikr,ikr)=Kd12(ikr,ikr)/2.d0
  enddo
  
  Ef_min=minval(Ef)
  
  ! calculate energy of A1 excitons with spin s=0 and s=1 *************************************************************
  Kernel=Ke11-Kd11+Kd12
  call evchf(Kernel,Ex_A1,Psi_A1)
  error=epihf(nkr,Kernel,Ex_A1,Psi_A1)
  if (error .ge. 1) stop "*** poor performance in calculating kernel eigenvalues ***"
   
  ! calculate energy of A2 excitons with spin s=0 *********************************************************************
  Kernel=Ke11+4.d0*Kx11-Kd11-Kd12
  call evchf(Kernel,Ex0_A2,Psi0_A2)
  error=epihf(nkr,Kernel,Ex0_A2,Psi0_A2)
  if (error .ge. 1) stop "*** poor performance in calculating kernel eigenvalues ***"
  
  ! calculate energy of A2 excitons with spin s=1 *********************************************************************
  Kernel=Ke11-Kd11-Kd12
  call evchf(Kernel,Ex1_A2,Psi1_A2)
  error=epihf(nkr,Kernel,Ex1_A2,Psi1_A2)
  if (error .ge. 1) stop "*** poor performance in calculating kernel eigenvalues ***"
    
  ! sort the subbands
  do ikr=ikr_high,ikr_low+1,-1
    do ikpr=ikr-1,ikr_low,-1
      
      if (Ex_A1(ikr) .lt. Ex_A1(ikpr)) then
        tmpr=Ex_A1(ikr)
        Psi_tmp=Psi_A1(:,ikr)
        Ex_A1(ikr)=Ex_A1(ikpr)
        Psi_A1(:,ikr)=Psi_A1(:,ikpr)
        Ex_A1(ikpr)=tmpr
        Psi_A1(:,ikpr)=Psi_tmp
      end if     
        
      if (Ex0_A2(ikr) .lt. Ex0_A2(ikpr)) then
        tmpr=Ex0_A2(ikr)
        Psi_tmp=Psi0_A2(:,ikr)
        Ex0_A2(ikr)=Ex0_A2(ikpr)
        Psi0_A2(:,ikr)=Psi0_A2(:,ikpr)
        Ex0_A2(ikpr)=tmpr
        Psi0_A2(:,ikpr)=Psi_tmp
      end if 
      
      if (Ex1_A2(ikr) .lt. Ex1_A2(ikpr)) then
        tmpr=Ex1_A2(ikr)
        Psi_tmp=Psi1_A2(:,ikr)
        Ex1_A2(ikr)=Ex1_A2(ikpr)
        Psi1_A2(:,ikr)=Psi1_A2(:,ikpr)
        Ex1_A2(ikpr)=tmpr
        Psi1_A2(:,ikpr)=Psi_tmp
      end if
    end do
  end do
  
  return
end