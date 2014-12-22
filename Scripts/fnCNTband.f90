subroutine fnCNTband()
  use comparams
  implicit none
  
  integer :: nkc, imin_sub
  integer :: i,j,mu,ik,tmpi
  integer, dimension(:), allocatable :: min_loc
  real*8 :: tmpr
  real*8, dimension(2) :: k,E1_tmp,E2_tmp
  real*8, dimension(:), allocatable :: k_vec,min_energy
  real*8, dimension(:,:,:), allocatable :: E_k
  complex*16, dimension(:,:,:), allocatable :: Cc_k,Cv_k
  complex*16, dimension(2) :: Cc_tmp,Cv_tmp
  
  
  ! calculate CNT energy dispersion.***********************************************************************************
  ikc_max=floor(pi/norm2(t_vec)/dk)
  ikc_min=-ikc_max
  nkc=2*ikc_max+1
  
  allocate(k_vec(ikc_min:ikc_max))
  allocate(E_k(1-Nu/2:Nu/2,ikc_min:ikc_max,2))
  allocate(Cc_k(1-Nu/2:Nu/2,ikc_min:ikc_max,2))
  allocate(Cv_k(1-Nu/2:Nu/2,ikc_min:ikc_max,2))
  allocate(min_loc(0:Nu/2))
  
  do ik=ikc_min,ikc_max
    k_vec(ik)=dble(ik)*dk
  end do
  
  do mu=1-Nu/2,Nu/2
    do ik=ikc_min,ikc_max
      k=dble(mu)*K1+dble(ik)*dk*K2
      call fnGrapheneEnergy(E_k(mu,ik,:),Cc_k(mu,ik,:),Cv_k(mu,ik,:),k)
    enddo
  enddo
  
  ! save the CNT energy dispersion*************************************************************************************
  do ik=ikc_min,ikc_max
     write(fh4,10, advance='no') k_vec(ik) 
     write(fh5,10, advance='no') k_vec(ik)
  end do
  write(fh4,10)
  write(fh5,10)
  
  do mu=1-Nu/2,Nu/2
    do ik=ikc_min,ikc_max
      write(fh4,10, advance='no') E_k(mu,ik,1) 
      write(fh5,10, advance='no') E_k(mu,ik,2) 
    end do
    write(fh4,10)
    write(fh5,10)
  enddo
  
  ! find the subbands with a minimum energy.***************************************************************************
  min_loc=minloc(E_k(0:Nu/2,:,1),2)
  imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
  allocate(min_sub(imin_sub))
  allocate(min_energy(imin_sub))
  
  i=1
  do mu=0,Nu/2
    if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
       min_sub(i)=mu
       min_energy(i)=minval(E_k(mu,:,1))
       i=i+1
    end if
  end do
  
  ! sort the subbands
  do i=imin_sub,2,-1
    do j=i-1,1,-1
      if (min_energy(i) .lt. min_energy(j)) then
        tmpr=min_energy(i)
        tmpi=min_sub(i)
        min_energy(i)=min_energy(j)
        min_sub(i)=min_sub(j)
        min_energy(j)=tmpr
        min_sub(j)=tmpi
      end if    
    end do
  end do
  
  ! find the max k-index that energy is below threshold energy (E_th).
  ik=0
  E1_tmp=(/ min_energy(i_sub),0 /)
  E2_tmp=(/ min_energy(i_sub),0 /)
  do while ((min(E1_tmp(1),E2_tmp(1))-min_energy(i_sub)) .le. E_th )
      k=dble(min_sub(i_sub))*K1+dble(ik)*dk*K2
      call fnGrapheneEnergy(E1_tmp,Cc_tmp,Cv_tmp,k)
      k=dble(min_sub(i_sub))*K1-dble(ik)*dk*K2
      call fnGrapheneEnergy(E2_tmp(:),Cc_tmp,Cv_tmp,k)
      ik=ik+1
  enddo
  
  ! set the index boundaries for some arrays and kernels. *************************************************************
  ik_max=ik                         !the higher limit of k-vector that is below E_th
  ik_min=-ik                        !the lower limit of k-vector that is below E_th
  iKcm_max=floor(Kcm_max/dk)        !the higher limit of center of mass wave vector that we calculate
  iKcm_min=-iKcm_max                !the lower limit of center of mass wave vector that we calculate
  ikr_high=iKcm_max-ik_min          !the maximum index that the relative wavenumber in the entire simulation.
  ikr_low=-ikr_high                 !the minimum index that the relative wavenumber in the entire simulation.
  ik_high=ikr_high+iKcm_max         !the maximum index that the wavenumber in the entire simulation.
  ik_low=-ik_high                   !the minimum index that the wavenumber in the entire simulation.
  iq_max=2*ikr_high                 !the higher limit of the index in v_FT and esp_q
  iq_min=-iq_max                    !the lower limit of the index in v_FT and esp_q
  
  ! calculate and save dispersion of the considered subbands **********************************************************
  deallocate(E_k)
  deallocate(Cc_k)
  deallocate(Cv_k)
  deallocate(k_vec)
  
  allocate(k_vec(ik_low:ik_high))
  allocate(E_k(2,ik_low:ik_high,2))
  allocate(Cc_k(2,ik_low:ik_high,2))
  allocate(Cv_k(2,ik_low:ik_high,2))
  
  do ik=ik_low,ik_high
    k_vec(ik)=dble(ik)*dk
  end do
  
  do ik=ik_low,ik_high
     write(fh6,10, advance='no') k_vec(ik) 
     write(fh7,10, advance='no') k_vec(ik)
  end do
  
  write(fh6,10)
  write(fh7,10)
  
  mu=min_sub(i_sub)
  do ik=ik_low,ik_high
    k=dble(mu)*K1+dble(ik)*dk*K2
    call fnGrapheneEnergy(E_k(1,ik,:),Cc_k(1,ik,:),Cv_k(1,ik,:),k)
    write(fh6,10, advance='no') E_k(1,ik,1) 
    write(fh7,10, advance='no') E_k(1,ik,2) 
  enddo
  
  write(fh6,10)
  write(fh7,10)
  
  mu=-min_sub(i_sub)
  do ik=ik_low,ik_high
    k=dble(mu)*K1+dble(ik)*dk*K2
    call fnGrapheneEnergy(E_k(2,ik,:),Cc_k(2,ik,:),Cv_k(2,ik,:),k)
    write(fh6,10, advance='no') E_k(2,ik,1) 
    write(fh7,10, advance='no') E_k(2,ik,2) 
  enddo
  
  
10 FORMAT (E16.8)
  
  continue
  return
end