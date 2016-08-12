module cnt_band_structure_mod
	implicit none
	private
	public :: cnt_band_structure, grapheneEnergy

contains
	subroutine cnt_band_structure()
		use comparams, only: currcnt
		use math_functions_mod, only: my_norm2
		use physical_constant_mod, only: pi
		use write_log_mod, only: writeLog
		implicit none
		
		integer :: nkc, imin_sub
		integer :: i,j,mu,ik,tmpi
		integer, dimension(:), allocatable :: min_loc
		real*8 :: tmpr
		real*8, dimension(2) :: k,E_tmp,E1_tmp,E2_tmp
		real*8, dimension(:), allocatable :: k_vec,min_energy
		real*8, dimension(:,:,:), allocatable :: E_k
		complex*16, dimension(:,:,:), allocatable :: Cc_k,Cv_k
		complex*16, dimension(2) :: Cc_tmp,Cv_tmp
		character(len=200) :: logInput
		
		
		! calculate CNT energy dispersion.***********************************************************************************
		currcnt%ikc_max=floor(pi/my_norm2(currcnt%t_vec)/currcnt%dk)
		currcnt%ikc_min=-currcnt%ikc_max
		nkc=2*currcnt%ikc_max+1
		
		allocate(k_vec(currcnt%ikc_min:currcnt%ikc_max))
		allocate(E_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(Cc_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(Cv_k(1-currcnt%Nu/2:currcnt%Nu/2,currcnt%ikc_min:currcnt%ikc_max,2))
		allocate(min_loc(0:currcnt%Nu/2))
		
		do ik=currcnt%ikc_min,currcnt%ikc_max
			k_vec(ik)=dble(ik)*currcnt%dk
		end do
		
		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
				call grapheneEnergy(currCNT,E_tmp,Cc_tmp,Cv_tmp,k)
				E_k(mu,ik,:) = E_tmp
				Cc_k(mu,ik,:) = Cc_tmp
				Cv_k(mu,ik,:) = Cv_tmp
			enddo
		enddo
		
		! save the CNT energy dispersion*************************************************************************************
		open(unit=100,file='CondBand.dat',status="unknown")
		open(unit=101,file='ValeBand.dat',status="unknown")
		
		do ik=currcnt%ikc_min,currcnt%ikc_max
			 write(100,'(E16.8)', advance='no') k_vec(ik) 
			 write(101,'(E16.8)', advance='no') k_vec(ik)
		end do
		
		write(100,*)
		write(101,*)
		
		do mu=1-currcnt%Nu/2,currcnt%Nu/2
			do ik=currcnt%ikc_min,currcnt%ikc_max
				write(100,'(E16.8)', advance='no') E_k(mu,ik,1) 
				write(101,'(E16.8)', advance='no') E_k(mu,ik,2) 
			end do
			write(100,*)
			write(101,*)
		enddo

		close(100)
		close(101)
		
		! find the subbands with a minimum energy.***************************************************************************
		min_loc=minloc(E_k(0:currcnt%Nu/2,:,1),2)
		imin_sub=count((min_loc .lt. nkc) .and. (min_loc .gt. 1))
		allocate(currcnt%min_sub(imin_sub))
		allocate(min_energy(imin_sub))
		
		i=1
		do mu=0,currcnt%Nu/2
			if ((min_loc(mu) .gt. 1) .and. (min_loc(mu) .lt. nkc)) then
			 currcnt%min_sub(i)=mu
			 min_energy(i)=minval(E_k(mu,:,1))
			 i=i+1
			end if
		end do
		
		! sort the subbands
		do i=imin_sub,2,-1
			do j=i-1,1,-1
			if (min_energy(i) .lt. min_energy(j)) then
				tmpr=min_energy(i)
				tmpi=currcnt%min_sub(i)
				min_energy(i)=min_energy(j)
				currcnt%min_sub(i)=currcnt%min_sub(j)
				min_energy(j)=tmpr
				currcnt%min_sub(j)=tmpi
			end if
			end do
		end do
		
		! find the max k-index that energy is below threshold energy (E_th).
		ik=0
		E1_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
		E2_tmp=(/ min_energy(currcnt%i_sub),0.d0 /)
		do while ((min(E1_tmp(1),E2_tmp(1))-min_energy(currcnt%i_sub)) .le. currcnt%E_th )
			k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E1_tmp,Cc_tmp,Cv_tmp,k)
			k=dble(currcnt%min_sub(currcnt%i_sub))*currcnt%K1-dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E2_tmp,Cc_tmp,Cv_tmp,k)
			ik=ik+1
		enddo
		
		! set the index boundaries for some arrays and kernels. *************************************************************
		currcnt%ik_max=ik											!the higher limit of k-vector that is below E_th
		currcnt%ik_min=-ik											!the lower limit of k-vector that is below E_th
		currcnt%iKcm_max=floor(currcnt%Kcm_max/currcnt%dk)			!the higher limit of center of mass wave vector that we calculate
		currcnt%iKcm_min=-currcnt%iKcm_max							!the lower limit of center of mass wave vector that we calculate
		currcnt%ikr_high=currcnt%iKcm_max-currcnt%ik_min			!the maximum index that the relative wavenumber in the entire simulation.
		currcnt%ikr_low=-currcnt%ikr_high							!the minimum index that the relative wavenumber in the entire simulation.
		currcnt%ik_high=currcnt%ikr_high+currcnt%iKcm_max			!the maximum index that the wavenumber in the entire simulation.
		currcnt%ik_low=-currcnt%ik_high								!the minimum index that the wavenumber in the entire simulation.
		currcnt%iq_max=max(2*currcnt%ikr_high,currcnt%ikc_max)		!the higher limit of the index in v_FT and esp_q
		currcnt%iq_min=-currcnt%iq_max								!the lower limit of the index in v_FT and esp_q

		! save the index boundaries to the log file. ************************************************************************
		call writeLog(new_line('A')//"Index boundaries *************************************")
		write(logInput,'(A, I0)') "ik_max=",currcnt%ik_max
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "iKcm_max=",currcnt%iKcm_max
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "ikr_high=",currcnt%ikr_high
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "ik_high=",currcnt%ik_high
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "iq_max=",currcnt%iq_max
		call writeLog(trim(logInput))
		
		! clean some variables that will be used for other purposes.
		deallocate(E_k)
		deallocate(Cc_k)
		deallocate(Cv_k)
		deallocate(k_vec)

		! calculate and save dispersion of the target subbands **********************************************************
		
		allocate(k_vec(currcnt%ik_low:currcnt%ik_high))
		allocate(E_k(2,currcnt%ik_low:currcnt%ik_high,2))
		allocate(Cc_k(2,currcnt%ik_low:currcnt%ik_high,2))
		allocate(Cv_k(2,currcnt%ik_low:currcnt%ik_high,2))
		
		do ik=currcnt%ik_low,currcnt%ik_high
			k_vec(ik)=dble(ik)*currcnt%dk
		end do
		
		open(unit=100,file='CondBand_Sub.dat',status="unknown")
		open(unit=101,file='ValeBand_Sub.dat',status="unknown")

		do ik=currcnt%ik_low,currcnt%ik_high
			write(100,'(E16.8)', advance='no') k_vec(ik) 
			write(101,'(E16.8)', advance='no') k_vec(ik)
		end do
		
		write(100,*)
		write(101,*)
		
		mu=currcnt%min_sub(currcnt%i_sub)
		do ik=currcnt%ik_low,currcnt%ik_high
			k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
			write(100,'(E16.8)', advance='no') E_tmp(1) 
			write(101,'(E16.8)', advance='no') E_tmp(2) 
		enddo
		
		write(100,*)
		write(101,*)
		
		mu=-currcnt%min_sub(currcnt%i_sub)
		do ik=currcnt%ik_low,currcnt%ik_high
			k=dble(mu)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
			write(100,'(E16.8)', advance='no') E_tmp(1) 
			write(101,'(E16.8)', advance='no') E_tmp(2) 
		enddo

		close(100)
		close(101)

		! write the index of the considered subband to the log file
		write(logInput,*) "min_sub(i_sub)=",currcnt%min_sub(currcnt%i_sub)
		call writeLog(trim(logInput))
		
		return
	end subroutine cnt_band_structure

	!**************************************************************************************************************************
	! subroutine to calculate Bloch functions and energy in graphene
	!**************************************************************************************************************************

	subroutine grapheneEnergy(currCNT,E,Cc,Cv,k)
		use physical_constant_mod, only: i1, t0
		use cnt_class, only: cnt

		type(cnt), intent(in) :: currCNT
		complex*16 :: f_k
		real*8, dimension(2), intent(in) :: k
		real*8, dimension(2), intent(out) :: E
		complex*16, dimension(2), intent(out) :: Cv
		complex*16, dimension(2), intent(out) :: Cc
  
		f_k=exp(i1*dcmplx(dot_product(k,(currCNT%a1+currCNT%a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(currCNT%a1-2.d0*currCNT%a2)/3.d0)))+exp(i1*dcmplx(dot_product(k,(currCNT%a2-2.d0*currCNT%a1)/3.d0)))
  
		E(1)=+t0*abs(f_k)
		E(2)=-t0*abs(f_k)
  
		Cc(1)=dcmplx(+1.d0/sqrt(2.d0))
		Cc(2)=dcmplx(+1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
		Cv(1)=dcmplx(+1.d0/sqrt(2.d0))
		Cv(2)=dcmplx(-1.d0/sqrt(2.d0)/abs(f_k))*conjg(f_k)
	end subroutine grapheneEnergy

end module cnt_band_structure_mod