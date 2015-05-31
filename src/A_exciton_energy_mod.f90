module A_exciton_energy_mod
	implicit none
	private
	public :: calculate_A_exciton_dispersion

	real*8, dimension(:), allocatable :: Ef
	complex*16, dimension(:), allocatable :: Psi_tmp
	complex*16, dimension(:,:), allocatable :: Kd11, Kd12, Kx11, Ke11, Kernel, Ktmp
	
	real*8, dimension(:), allocatable :: Ex_A1, Ex0_A2, Ex1_A2
	complex*16, dimension(:,:), allocatable :: Psi_A1, Psi0_A2, Psi1_A2 !the first index is ikr, the scond index is the subband

contains
	subroutine calculate_A_exciton_dispersion()
		use comparams, only: currcnt
		use write_log_mod, only: writeLog
		
		integer :: iKcm, ikr, ikpr
		real*8 :: Ef_min
		character(len=200) :: logInput

		allocate(currcnt%Psi_A1(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi0_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi1_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex_A1(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex0_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex1_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))


		! find the number of exciton bands below the free electron level
		currcnt%nX=0
		iKcm=0
		Ef_min = 0.d0
		call calculate_exciton_energy(Ef_min, iKcm)

		do ikr=currcnt%ikr_low,currcnt%ikr_high
			if (Ef_min .ge. Ex0_A2(ikr)) then
				currcnt%nX=currcnt%nX+1
			endif
		enddo

		write(logInput,*) "nX=",currcnt%nX
		call writeLog(new_line('A')//trim(logInput))
		
		call writeLog(new_line('A')//"Calculating exciton dispersion ********************************")
		

		open(unit=100,file='Ex_A1.dat',status="unknown")
		open(unit=101,file='Ex0_A2.dat',status="unknown")
		open(unit=102,file='Ex1_A2.dat',status="unknown")
		open(unit=103,file='Psi_A1.dat',status="unknown")
		open(unit=104,file='Psi0_A2.dat',status="unknown")
		open(unit=105,file='Psi1_A2.dat',status="unknown")

		do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
			write(logInput,*) "iKcm=", iKcm
			call writeLog(trim(logInput))
			
			call calculate_exciton_energy(Ef_min, iKcm)

			currcnt%Ex_A1(:,iKcm) = Ex_A1
			currcnt%Ex0_A2(:,iKcm) = Ex0_A2
			currcnt%Ex1_A2(:,iKcm) = Ex1_A2
			currcnt%Psi_A1(:,:,iKcm) = Psi_A1
			currcnt%Psi0_A2(:,:,iKcm) = Psi0_A2
			currcnt%Psi1_A2(:,:,iKcm) = Psi1_A2
			
			! save exciton energy and wavefunction
! 			do ikr=(currcnt%ikr_low),(currcnt%ikr_low+currcnt%nX-1)
			do ikr=(currcnt%ikr_low),(currcnt%ikr_high)
				write(100,'(E16.8)', advance='no') Ex_A1(ikr)
				write(101,'(E16.8)', advance='no') Ex0_A2(ikr)
				write(102,'(E16.8)', advance='no') Ex1_A2(ikr)
				do ikpr=currcnt%ikr_low,currcnt%ikr_high
					write(103,'(E16.8,E16.8)', advance='no') Psi_A1(ikpr,ikr)
					write(104,'(E16.8,E16.8)', advance='no') Psi0_A2(ikpr,ikr)
					write(105,'(E16.8,E16.8)', advance='no') Psi1_A2(ikpr,ikr)
				enddo
			enddo
			write(100,*)
			write(101,*)
			write(102,*)
			write(103,*)
			write(104,*)
			write(105,*)
		
		enddo
		
		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		close(105)
		  
		return
	end subroutine calculate_A_exciton_dispersion


	!**************************************************************************************************************************
	! subroutine to calculate A-type exciton energy with center of mass k vector iKcm
	!**************************************************************************************************************************

	subroutine calculate_exciton_energy(Ef_min, iKcm)
		use comparams, only: currcnt
		use math_functions_mod, only: eig
		
		integer :: nkr
		integer :: iKcm,mu_cm,mu_kr
		integer :: ikr, ikc, ikv, ikpr, ikpc, ikpv
		real*8 :: tmpr, Ef_min
		real*8, dimension(2) :: Kcm
		

		mu_cm=0
		mu_kr=currcnt%min_sub(currcnt%i_sub)
		nkr=currcnt%ikr_high-currcnt%ikr_low+1

		if (.not. allocated(Psi_A1)) allocate(Psi_A1(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if (.not. allocated(Psi0_A2)) allocate(Psi0_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if (.not. allocated(Psi1_A2)) allocate(Psi1_A2(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if (.not. allocated(Ex_A1)) allocate(Ex_A1(currcnt%ikr_low:currcnt%ikr_high))
		if (.not. allocated(Ex0_A2)) allocate(Ex0_A2(currcnt%ikr_low:currcnt%ikr_high))
		if (.not. allocated(Ex1_A2)) allocate(Ex1_A2(currcnt%ikr_low:currcnt%ikr_high))
		
		if(.not. allocated(Kd11)) allocate(Kd11(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Kd12)) allocate(Kd12(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Kx11)) allocate(Kx11(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ke11)) allocate(Ke11(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Kernel)) allocate(Kernel(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ktmp)) allocate(Ktmp(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ef)) allocate(Ef(currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Psi_tmp)) allocate(Psi_tmp(currcnt%ikr_low:currcnt%ikr_high))
		
		Kcm=dble(mu_cm)*currcnt%K1+dble(iKcm)*currcnt%dk*currcnt%K2
			
		Ke11=0.d0*Ke11
		Kd11=0.d0*Kd11
		Kd12=0.d0*Kd12
		Kx11=0.d0*Kx11
		
		! calculate the kernel matrices
		do ikr=currcnt%ikr_low,currcnt%ikr_high
			ikc=ikr+iKcm
			ikv=ikr-iKcm
			Ke11(ikr,ikr)=dcmplx(currcnt%Ek(1,ikc,1)-currcnt%Ek(1,ikv,2)+currcnt%Sk(1,ikc,1)-currcnt%Sk(1,ikv,2))
			Ef(ikr)=currcnt%Ek(1,ikc,1)-currcnt%Ek(1,ikv,2)+currcnt%Sk(1,ikc,1)-currcnt%Sk(1,ikv,2)
			do ikpr=ikr,currcnt%ikr_high
				ikpc=ikpr+iKcm
				ikpv=ikpr-iKcm
				Kd11(ikr,ikpr)=(conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(1,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(1,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(1,ikpv,2))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(1,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(1,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(1,ikpv,2)))/dcmplx(currcnt%kappa*currcnt%eps_q(0,ikr-ikpr))
					
				Kx11(ikr,ikpr)=(conjg(currcnt%Cc(1,ikc,1))*currcnt%Cv(1,ikv,1)*currcnt%v_FT(0,2*iKcm,1,1)*currcnt%Cc(1,ikpc,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,1))*currcnt%Cv(1,ikv,1)*currcnt%v_FT(0,2*iKcm,1,2)*currcnt%Cc(1,ikpc,2)*conjg(currcnt%Cv(1,ikpv,2))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cv(1,ikv,2)*currcnt%v_FT(0,2*iKcm,2,1)*currcnt%Cc(1,ikpc,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cv(1,ikv,2)*currcnt%v_FT(0,2*iKcm,2,2)*currcnt%Cc(1,ikpc,2)*conjg(currcnt%Cv(1,ikpv,2)))
			enddo
			
			do ikpr=currcnt%ikr_high,-ikr,-1
				ikpc=ikpr+iKcm
				ikpv=ikpr-iKcm
				Kd12(ikr,-ikpr)=(conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(2,ikpc,1)*currcnt%v_FT(2*mu_kr,ikr-ikpr,1,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								 conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(2,ikpc,1)*currcnt%v_FT(2*mu_kr,ikr-ikpr,1,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(2,ikpv,2))+ &
								 conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(2,ikpc,2)*currcnt%v_FT(2*mu_kr,ikr-ikpr,2,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								 conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(2,ikpc,2)*currcnt%v_FT(2*mu_kr,ikr-ikpr,2,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(2,ikpv,2)))/dcmplx(currcnt%kappa*currcnt%eps_q(2*mu_kr,ikr-ikpr))
			enddo
		enddo
		
		! when running the code in release mode there the next few lines generated a stack overflow error which forced me to use a dummy variable to resolve this issue.
		Ktmp=conjg(transpose(Kd11))
		Kd11=Kd11+Ktmp
		Ktmp=conjg(transpose(Kx11))
		Kx11=Kx11+Ktmp
		Ktmp=conjg(transpose(Kd12))
		Kd12=Kd12+Ktmp
		
		do ikr=currcnt%ikr_low,currcnt%ikr_high
			Kd11(ikr,ikr)=Kd11(ikr,ikr)/2.d0
			Kx11(ikr,ikr)=Kx11(ikr,ikr)/2.d0
			Kd12(ikr,ikr)=Kd12(ikr,ikr)/2.d0
		enddo
		
		Ef_min=minval(Ef)
		
		! calculate energy of A1 excitons with spin s=0 and s=1 *************************************************************
		Kernel=Ke11-Kd11+Kd12
		call eig(nkr,Kernel,Psi_A1,Ex_A1)
		 
		! calculate energy of A2 excitons with spin s=0 *********************************************************************
		Kernel=Ke11+4.d0*Kx11-Kd11-Kd12
		call eig(nkr,Kernel,Psi0_A2,Ex0_A2)
		
		! calculate energy of A2 excitons with spin s=1 *********************************************************************
		Kernel=Ke11-Kd11-Kd12
		call eig(nkr,Kernel,Psi1_A2,Ex1_A2)
			
		!sort the subbands
		do ikr=currcnt%ikr_high,currcnt%ikr_low+1,-1
			do ikpr=ikr-1,currcnt%ikr_low,-1
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
	end subroutine calculate_exciton_energy

end module A_exciton_energy_mod