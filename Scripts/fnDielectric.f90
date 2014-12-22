subroutine fnDielectric()
  use comparams
  implicit none
  
  integer :: ik, mu_k, iq, mu_q, ikq, mu_kq
  integer :: tmpi,i,u,s
  real*8, dimension(2) :: deltaR
  real*8, dimension(2) :: k, q, kq
  real*8, dimension(:,:), allocatable :: PI_q,v_q
  real*8, dimension(:,:,:), allocatable :: E_k
  real*8, dimension(:,:,:), allocatable :: q_vec
  !real*8, dimension(:), allocatable :: q_vec
  real*8, dimension(2) :: Ek_tmp, Ekq_tmp
  complex*16, dimension(2) :: Cc_tmp, Cv_tmp, Ccq_tmp, Cvq_tmp
  complex*16, dimension(:,:,:), allocatable :: Cc_k, Cv_k
  
  
  ! initialize variables.**********************************************************************************************
  allocate(PI_q(1-Nu:Nu-1,iq_min:iq_max))
  allocate(v_q(1-Nu:Nu-1,iq_min:iq_max))
  allocate(eps_q(1-Nu:Nu-1,iq_min:iq_max))
  allocate(v_FT(1-Nu:Nu-1,iq_min:iq_max,2,2))
  allocate(q_vec(1-Nu:Nu-1,iq_min:iq_max,2))
  !allocate(q_vec(iq_min:iq_max))
  
  allocate(Cc_k(2-3*Nu/2:3*Nu/2-1,iq_min+ikc_min:iq_max+ikc_max,2))
  allocate(Cv_k(2-3*Nu/2:3*Nu/2-1,iq_min+ikc_min:iq_max+ikc_max,2))
  allocate(E_k(2-3*Nu/2:3*Nu/2-1,iq_min+ikc_min:iq_max+ikc_max,2))
  
  do mu_q=1-Nu,Nu-1
    do iq=iq_min,iq_max
      PI_q(mu_q,iq)=(0.d0,0.d0)
      v_q(mu_q,iq)=(0.d0,0.d0)
      eps_q(mu_q,iq)=(0.d0,0.d0)
      v_FT(mu_q,iq,1,1)=(0.d0,0.d0)
      v_FT(mu_q,iq,1,2)=(0.d0,0.d0)
      v_FT(mu_q,iq,2,1)=(0.d0,0.d0)
      v_FT(mu_q,iq,2,2)=(0.d0,0.d0)
    enddo
  enddo
  
  if (flg_dielectric .eq. .true.) then
  ! calculate PI(q).***************************************************************************************************
    do mu_kq=2-3*Nu/2,3*Nu/2-1
      do ikq=(iq_min+ikc_min),(iq_max+ikc_max)
        kq=dble(mu_kq)*K1+dble(ikq)*dk*K2
        call fnGrapheneEnergy(E_k(mu_kq,ikq,:),Cc_k(mu_kq,ikq,:),Cv_k(mu_kq,ikq,:),kq)
      enddo
    enddo
    
    do mu_q=1-Nu,Nu-1
      print *,"mu_q=",mu_q
      do iq=iq_min,iq_max
        do mu_k=1-Nu/2,Nu/2
          do ik=ikc_min,ikc_max
            ikq=ik+iq
            mu_kq=mu_k+mu_q
            PI_q(mu_q,iq)=PI_q(mu_q,iq)+ &
                (abs(dot_product(Cv_k(mu_k,ik,:),Cc_k(mu_kq,ikq,:))))**2/(E_k(mu_kq,ikq,1)-E_k(mu_k,ik,2))+ &
                (abs(dot_product(Cc_k(mu_k,ik,:),Cv_k(mu_kq,ikq,:))))**2/(E_k(mu_k,ik,1)-E_k(mu_kq,ikq,2))
          enddo
        enddo
      enddo
    enddo
  
  ! calculate v_FT and eps_q.******************************************************************************************
    do mu_q=1-Nu,Nu-1
      do iq=iq_min,iq_max
        q_vec(mu_q,iq,1)=dble(mu_q)*K1(1)+dble(iq)*dk*K2(1)
        q_vec(mu_q,iq,2)=dble(mu_q)*K1(2)+dble(iq)*dk*K2(2)
      enddo
    enddo
    
    do u=-nr,nr
      print *,"u=",u
      do s=1,Nu 
        deltaR=dble(u)*t_vec+posAA(s,:)
        v_FT(:,:,1,1)=v_FT(:,:,1,1)+exp(i1*dcmplx(q_vec(:,:,1)*deltaR(1)+q_vec(:,:,2)*deltaR(2)))* &
                dcmplx(Upp/sqrt((4.d0*pi*eps0/q0**2*Upp*norm2(deltaR))**2+1))
          
        deltaR=dble(u)*t_vec+posBA(s,:)
        v_FT(:,:,1,2)=v_FT(:,:,1,2)+exp(i1*dcmplx(q_vec(:,:,1)*deltaR(1)+q_vec(:,:,2)*deltaR(2)))* &
                dcmplx(Upp/sqrt((4.d0*pi*eps0/q0**2*Upp*norm2(deltaR))**2+1))
      enddo
    enddo
    
    v_FT(:,:,2,2)=v_FT(:,:,1,1)
    
    do mu_q=1-Nu,Nu-1
      do iq=iq_min,iq_max
        v_FT(mu_q,iq,2,1)=v_FT(-mu_q,-iq,1,2)
      enddo
    enddo
    
    do mu_q=1-Nu,Nu-1
      do iq=iq_min,iq_max
        v_q(mu_q,iq)=(1.d0/4.d0)*dreal(sum(v_FT(mu_q,iq,:,:)))
      enddo
    enddo
    
    v_FT=v_FT/dcmplx((ikc_max-ikc_min+1)*Nu)
    v_q=v_q/dble((ikc_max-ikc_min+1)*Nu)
    
    eps_q=1.d0+v_q*PI_q
  
  ! calculate PI(q).***************************************************************************************************    
    !do ik=ikc_min,ikc_max
    !  print *, "ik=", ik
    !  do mu_k=1-Nu/2,Nu/2
    !    k=dble(mu_k)*K1+dble(ik)*dk*K2
    !    call fnGrapheneEnergy(Ek_tmp,Cc_tmp,Cv_tmp,k)
    !    do iq=iq_min,iq_max
    !      do mu_q=1-Nu,Nu-1
    !        mu_kq=mu_k+mu_q
    !        kq=dble(mu_kq)*K1+dble(ik+iq)*dk*K2
    !        call fnGrapheneEnergy(Ekq_tmp,Ccq_tmp,Cvq_tmp,kq)
    !        PI_q(mu_q,iq)=PI_q(mu_q,iq)+ &
    !            ((abs(dot_product(Cv_tmp,Ccq_tmp)))**2)/(Ekq_tmp(1)-Ek_tmp(2))+ &
    !            ((abs(dot_product(Cc_tmp,Cvq_tmp)))**2)/(Ek_tmp(1)-Ekq_tmp(2))
    !      enddo
    !    enddo
    !  enddo
    !enddo
  
  ! calculate v_FT and eps_q.******************************************************************************************
    !do iq=iq_min,iq_max
    !  q_vec(iq)=dble(iq)*dk
    !enddo
    !    
    !do mu_q=1-Nu,Nu-1
    !  print *,"mu_q=",mu_q
    !  do u=-nr,nr
    !    do s=1,Nu
    !      deltaR=dble(u)*t_vec+posAA(s,:)
    !      v_FT(mu_q,:,1,1)=v_FT(mu_q,:,1,1)+exp(i1*(q_vec*dot_product(K2,deltaR)+dble(mu_q)*dot_product(K1,deltaR)))* &
    !          dcmplx(Upp/sqrt((4.d0*pi*eps0/(q0**2)*Upp*norm2(deltaR))**2+1.d0))
    !      
    !      deltaR=dble(u)*t_vec+posBA(s,:)
    !      v_FT(mu_q,:,1,2)=v_FT(mu_q,:,1,2)+exp(i1*(q_vec*dot_product(K2,deltaR)+dble(mu_q)*dot_product(K1,deltaR)))* &
    !          dcmplx(Upp/sqrt((4.d0*pi*eps0/(q0**2)*Upp*norm2(deltaR))**2+1.d0))
    !      
    !      deltaR=dble(u)*t_vec+posAB(s,:)
    !      v_FT(mu_q,:,2,1)=v_FT(mu_q,:,2,1)+exp(i1*(q_vec*dot_product(K2,deltaR)+dble(mu_q)*dot_product(K1,deltaR)))* &
    !          dcmplx(Upp/sqrt((4.d0*pi*eps0/(q0**2)*Upp*norm2(deltaR))**2+1.d0))
    !      
    !      deltaR=dble(u)*t_vec+posBB(s,:)
    !      v_FT(mu_q,:,2,2)=v_FT(mu_q,:,2,2)+exp(i1*(q_vec*dot_product(K2,deltaR)+dble(mu_q)*dot_product(K1,deltaR)))* &
    !          dcmplx(Upp/sqrt((4.d0*pi*eps0/(q0**2)*Upp*norm2(deltaR))**2+1.d0))
    !    enddo
    !  enddo    
    !  v_q(mu_q,:)=1.d0/4.d0*dreal(v_FT(mu_q,:,1,1)+v_FT(mu_q,:,1,2)+v_FT(mu_q,:,2,1)+v_FT(mu_q,:,2,2))
    !enddo
    !
    !v_FT=v_FT/dcmplx((2*ikc_max+1)*Nu)
    !v_q=v_q/dble((2*ikc_max+1)*Nu)
    !
    !eps_q=1.d0+v_q*PI_q
      
      
  ! save eps_q, and v_FT.**********************************************************************************************
  
    do mu_q=1-Nu,Nu-1
      do iq=iq_min,iq_max
        write(fh8,10, advance='no') PI_q(mu_q,iq)
        write(fh9,11, advance='no') v_FT(mu_q,iq,1,1) 
        write(fh9,11, advance='no') v_FT(mu_q,iq,1,2) 
        write(fh9,11, advance='no') v_FT(mu_q,iq,2,1) 
        write(fh9,11, advance='no') v_FT(mu_q,iq,2,2)
        write(fh10,10, advance='no') eps_q(mu_q,iq) 
      end do
      
      write(fh8,10)
      write(fh9,10)
      write(fh10,10)
    enddo
  
  else
    do mu_q=1-Nu,Nu-1
      do iq=iq_min,iq_max
        read(fh8,10, advance='no') PI_q(mu_q,iq) 
        read(fh9,11, advance='no') v_FT(mu_q,iq,1,1) 
        read(fh9,11, advance='no') v_FT(mu_q,iq,1,2) 
        read(fh9,11, advance='no') v_FT(mu_q,iq,2,1) 
        read(fh9,11, advance='no') v_FT(mu_q,iq,2,2)
        read(fh10,10, advance='no') eps_q(mu_q,iq)
      end do
      
      read(fh8,10)
      read(fh9,10)
      read(fh10,10)
    enddo  
  endif
  
10 FORMAT (E16.8)  
11 FORMAT (E16.8,E16.8)  
   
  return
end