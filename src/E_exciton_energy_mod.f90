module E_exciton_energy_mod
	implicit none
	private
	public :: calculate_E_exciton_dispersion

	real*8, dimension(:), allocatable :: Ef
	complex*16, dimension(:), allocatable :: Psi_tmp
	complex*16, dimension(:,:), allocatable :: Kd, Kx, Ke, Kernel, Ktmp

	real*8, dimension(:), allocatable :: Ex0_Ep, Ex1_Ep, Ex0_Em, Ex1_Em
	complex*16, dimension(:,:), allocatable :: Psi0_Ep, Psi1_Ep, Psi0_Em, Psi1_Em !the first index is ikr, the scond index is the subband

contains
	subroutine calculate_E_exciton_dispersion()
		use comparams, only: currcnt
		use write_log_mod, only: writeLog
		
		integer :: iKcm, ikr, ikpr
		real*8 :: Ef_min
		character(len=200) :: logInput

		allocate(currcnt%Psi0_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi0_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi1_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Psi1_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex0_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex0_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex1_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))
		allocate(currcnt%Ex1_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%iKcm_min:currcnt%iKcm_max))

		! find the number of exciton bands below the free electron level
		currcnt%nX=0
		iKcm=0
		Ef_min = 0.d0
		call calculate_exciton_energy(Ef_min, iKcm)

		do ikr=currcnt%ikr_low,currcnt%ikr_high
			if (Ef_min .ge. Ex0_Ep(ikr)) then
				currcnt%nX=currcnt%nX+1
			endif
		enddo

		write(logInput,*) "nX=",currcnt%nX
		call writeLog(new_line('A')//trim(logInput))
		
		call writeLog(new_line('A')//"Calculating exciton dispersion ********************************")
		

		open(unit=100,file='Ex0_Ep.dat',status="unknown")
		open(unit=101,file='Ex0_Em.dat',status="unknown")
		open(unit=102,file='Ex1_Ep.dat',status="unknown")
		open(unit=103,file='Ex1_Em.dat',status="unknown")
		open(unit=104,file='Psi0_Ep.dat',status="unknown")
		open(unit=105,file='Psi0_Em.dat',status="unknown")
		open(unit=106,file='Psi1_Ep.dat',status="unknown")
		open(unit=107,file='Psi1_Em.dat',status="unknown")

		do iKcm=currcnt%iKcm_min,currcnt%iKcm_max
			write(logInput,*) "iKcm=", iKcm
			call writeLog(trim(logInput))
			
			call calculate_exciton_energy(Ef_min, iKcm)

			currcnt%Ex0_Ep(:,iKcm) = Ex0_Ep
			currcnt%Ex0_Em(:,iKcm) = Ex0_Em
			currcnt%Ex1_Ep(:,iKcm) = Ex1_Ep
			currcnt%Ex1_Em(:,iKcm) = Ex1_Em
			currcnt%Psi0_Ep(:,:,iKcm) = Psi0_Ep
			currcnt%Psi0_Em(:,:,iKcm) = Psi0_Em
			currcnt%Psi1_Ep(:,:,iKcm) = Psi1_Ep
			currcnt%Psi1_Em(:,:,iKcm) = Psi1_Em
			
			! save exciton energy and wavefunction
			do ikr=(currcnt%ikr_low),(currcnt%ikr_low+currcnt%nX-1)
				write(100,'(E16.8)', advance='no') Ex0_Ep(ikr)
				write(101,'(E16.8)', advance='no') Ex0_Em(ikr)
				write(102,'(E16.8)', advance='no') Ex1_Ep(ikr)
				write(103,'(E16.8)', advance='no') Ex1_Em(ikr)
				do ikpr=currcnt%ikr_low,currcnt%ikr_high
					write(104,'(E16.8,E16.8)', advance='no') Psi0_Ep(ikpr,ikr)
					write(105,'(E16.8,E16.8)', advance='no') Psi0_Em(ikpr,ikr)
					write(106,'(E16.8,E16.8)', advance='no') Psi1_Ep(ikpr,ikr)
					write(107,'(E16.8,E16.8)', advance='no') Psi1_Em(ikpr,ikr)
				enddo
			enddo
			write(100,*)
			write(101,*)
			write(102,*)
			write(103,*)
			write(104,*)
			write(105,*)
			write(106,*)
			write(107,*)
		
		enddo
		
		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		close(105)
		close(106)
		close(107)
		  
		return
	end subroutine calculate_E_exciton_dispersion

	!**************************************************************************************************************************
	! subroutine to calculate E-type exciton energy with center of mass k vector iKcm
	!**************************************************************************************************************************

	subroutine calculate_exciton_energy(Ef_min, iKcm)
		use comparams, only: currcnt
		use math_functions_mod, only: eig
		
		integer, intent(in) :: iKcm
		integer :: nkr
		integer :: mu_cm
		integer :: mu_c, mu_v
		integer :: ikr, ikc, ikv, ikpr, ikpc, ikpv
		real*8 :: tmpr, Ef_min
		
		nkr=currcnt%ikr_high-currcnt%ikr_low+1

		! E-type exciton
		if(.not. allocated(Psi0_Ep)) allocate(Psi0_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Psi0_Em)) allocate(Psi0_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Psi1_Ep)) allocate(Psi1_Ep(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Psi1_Em)) allocate(Psi1_Em(currcnt%ikr_low:currcnt%ikr_high, currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ex0_Ep)) allocate(Ex0_Ep(currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ex0_Em)) allocate(Ex0_Em(currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ex1_Ep)) allocate(Ex1_Ep(currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ex1_Em)) allocate(Ex1_Em(currcnt%ikr_low:currcnt%ikr_high))
		
		if(.not. allocated(Kd)) allocate(Kd(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Kx)) allocate(Kx(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ke)) allocate(Ke(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Kernel)) allocate(Kernel(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ktmp)) allocate(Ktmp(currcnt%ikr_low:currcnt%ikr_high,currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Ef)) allocate(Ef(currcnt%ikr_low:currcnt%ikr_high))
		if(.not. allocated(Psi_tmp)) allocate(Psi_tmp(currcnt%ikr_low:currcnt%ikr_high))
		
		

		! calculate the kernel matrices for E+ exciton **********************************************************************
		mu_c = +currcnt%min_sub(currcnt%i_sub)
		mu_v = -currcnt%min_sub(currcnt%i_sub)
		mu_cm = (mu_c-mu_v)/2

		Ke=0.d0*Ke
		Kd=0.d0*Kd
		Kx=0.d0*Kx

		do ikr=currcnt%ikr_low,currcnt%ikr_high
			ikc = ikr+iKcm
			ikv = ikr-iKcm

			Ke(ikr,ikr)=dcmplx(currcnt%Ek(1,ikc,1)-currcnt%Ek(2,ikv,2)+currcnt%Sk(1,ikc,1)-currcnt%Sk(2,ikv,2))
			Ef(ikr)=currcnt%Ek(1,ikc,1)-currcnt%Ek(2,ikv,2)+currcnt%Sk(1,ikc,1)-currcnt%Sk(2,ikv,2)
			do ikpr=ikr,currcnt%ikr_high
				ikpc=ikpr+iKcm
				ikpv=ikpr-iKcm
				Kd(ikr,ikpr)=(conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(1,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,1)*currcnt%Cv(2,ikv,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,1))*currcnt%Cc(1,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,2)*currcnt%Cv(2,ikv,2)*conjg(currcnt%Cv(2,ikpv,2))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(1,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,1)*currcnt%Cv(2,ikv,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cc(1,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,2)*currcnt%Cv(2,ikv,2)*conjg(currcnt%Cv(2,ikpv,2)))/dcmplx(currcnt%kappa*currcnt%eps_q(0,ikr-ikpr))
					
				Kx(ikr,ikpr)=(conjg(currcnt%Cc(1,ikc,1))*currcnt%Cv(2,ikv,1)*currcnt%v_FT(2*mu_cm,2*iKcm,1,1)*currcnt%Cc(1,ikpc,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,1))*currcnt%Cv(2,ikv,1)*currcnt%v_FT(2*mu_cm,2*iKcm,1,2)*currcnt%Cc(1,ikpc,2)*conjg(currcnt%Cv(2,ikpv,2))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cv(2,ikv,2)*currcnt%v_FT(2*mu_cm,2*iKcm,2,1)*currcnt%Cc(1,ikpc,1)*conjg(currcnt%Cv(2,ikpv,1))+ &
								conjg(currcnt%Cc(1,ikc,2))*currcnt%Cv(2,ikv,2)*currcnt%v_FT(2*mu_cm,2*iKcm,2,2)*currcnt%Cc(1,ikpc,2)*conjg(currcnt%Cv(2,ikpv,2)))
			enddo
		enddo
		
		! when running the code in release mode there the next few lines generated a stack overflow error which forced me to use a dummy variable to resolve this issue.
		Ktmp=conjg(transpose(Kd))
		Kd=Kd+Ktmp
		Ktmp=conjg(transpose(Kx))
		Kx=Kx+Ktmp
		
		do ikr=currcnt%ikr_low,currcnt%ikr_high
			Kd(ikr,ikr)=Kd(ikr,ikr)/2.d0
			Kx(ikr,ikr)=Kx(ikr,ikr)/2.d0
		enddo
		
		Ef_min=minval(Ef)
		
		! calculate energy of E+ excitons with spin s=0
		Kernel=Ke+dcmplx(2.d0)*Kx-Kd
		call eig(nkr,Kernel,Psi0_Ep,Ex0_Ep)
		 
		! calculate energy of E+ excitons with spin s=1
		Kernel=Ke-Kd
		call eig(nkr,Kernel,Psi1_Ep,Ex1_Ep)



		! calculate the kernel matrices for E- exciton *******************************************************************
		mu_c = -currcnt%min_sub(currcnt%i_sub)
		mu_v = +currcnt%min_sub(currcnt%i_sub)
		mu_cm = (mu_c-mu_v)/2

		Ke=0.d0*Ke
		Kd=0.d0*Kd
		Kx=0.d0*Kx

		do ikr=currcnt%ikr_low,currcnt%ikr_high
			ikc = ikr+iKcm
			ikv = ikr-iKcm

			Ke(ikr,ikr)=dcmplx(currcnt%Ek(2,ikc,1)-currcnt%Ek(1,ikv,2)+currcnt%Sk(2,ikc,1)-currcnt%Sk(1,ikv,2))
			Ef(ikr)=currcnt%Ek(2,ikc,1)-currcnt%Ek(1,ikv,2)+currcnt%Sk(2,ikc,1)-currcnt%Sk(1,ikv,2)
			do ikpr=ikr,currcnt%ikr_high
				ikpc=ikpr+iKcm
				ikpv=ikpr-iKcm
				Kd(ikr,ikpr)=(conjg(currcnt%Cc(2,ikc,1))*currcnt%Cc(2,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(2,ikc,1))*currcnt%Cc(2,ikpc,1)*currcnt%v_FT(0,ikr-ikpr,1,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(1,ikpv,2))+ &
								conjg(currcnt%Cc(2,ikc,2))*currcnt%Cc(2,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,1)*currcnt%Cv(1,ikv,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(2,ikc,2))*currcnt%Cc(2,ikpc,2)*currcnt%v_FT(0,ikr-ikpr,2,2)*currcnt%Cv(1,ikv,2)*conjg(currcnt%Cv(1,ikpv,2)))/dcmplx(currcnt%kappa*currcnt%eps_q(0,ikr-ikpr))
					
				Kx(ikr,ikpr)=(conjg(currcnt%Cc(2,ikc,1))*currcnt%Cv(1,ikv,1)*currcnt%v_FT(2*mu_cm,2*iKcm,1,1)*currcnt%Cc(2,ikpc,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(2,ikc,1))*currcnt%Cv(1,ikv,1)*currcnt%v_FT(2*mu_cm,2*iKcm,1,2)*currcnt%Cc(2,ikpc,2)*conjg(currcnt%Cv(1,ikpv,2))+ &
								conjg(currcnt%Cc(2,ikc,2))*currcnt%Cv(1,ikv,2)*currcnt%v_FT(2*mu_cm,2*iKcm,2,1)*currcnt%Cc(2,ikpc,1)*conjg(currcnt%Cv(1,ikpv,1))+ &
								conjg(currcnt%Cc(2,ikc,2))*currcnt%Cv(1,ikv,2)*currcnt%v_FT(2*mu_cm,2*iKcm,2,2)*currcnt%Cc(2,ikpc,2)*conjg(currcnt%Cv(1,ikpv,2)))
			enddo
		enddo
		
		! when running the code in release mode there the next few lines generated a stack overflow error which forced me to use a dummy variable to resolve this issue.
		Ktmp=conjg(transpose(Kd))
		Kd=Kd+Ktmp
		Ktmp=conjg(transpose(Kx))
		Kx=Kx+Ktmp
		
		do ikr=currcnt%ikr_low,currcnt%ikr_high
			Kd(ikr,ikr)=Kd(ikr,ikr)/2.d0
			Kx(ikr,ikr)=Kx(ikr,ikr)/2.d0
		enddo
		
		Ef_min=minval(Ef)
		
		! calculate energy of E- excitons with spin s=0
		Kernel=Ke+dcmplx(2.d0)*Kx-Kd
		call eig(nkr,Kernel,Psi0_Em,Ex0_Em)
		 
		! calculate energy of E- excitons with spin s=1
		Kernel=Ke-Kd
		call eig(nkr,Kernel,Psi1_Em,Ex1_Em)
			
		!sort the subbands ****************************************************************************************
		do ikr=currcnt%ikr_high,currcnt%ikr_low+1,-1
			do ikpr=ikr-1,currcnt%ikr_low,-1
				if (Ex0_Ep(ikr) .lt. Ex0_Ep(ikpr)) then
					tmpr=Ex0_Ep(ikr)
					Psi_tmp=Psi0_Ep(:,ikr)
					Ex0_Ep(ikr)=Ex0_Ep(ikpr)
					Psi0_Ep(:,ikr)=Psi0_Ep(:,ikpr)
					Ex0_Ep(ikpr)=tmpr
					Psi0_Ep(:,ikpr)=Psi_tmp
				end if 
				
				if (Ex1_Ep(ikr) .lt. Ex1_Ep(ikpr)) then
					tmpr=Ex1_Ep(ikr)
					Psi_tmp=Psi1_Ep(:,ikr)
					Ex1_Ep(ikr)=Ex1_Ep(ikpr)
					Psi1_Ep(:,ikr)=Psi1_Ep(:,ikpr)
					Ex1_Ep(ikpr)=tmpr
					Psi1_Ep(:,ikpr)=Psi_tmp
				end if 	
				
				if (Ex0_Em(ikr) .lt. Ex0_Em(ikpr)) then
					tmpr=Ex0_Em(ikr)
					Psi_tmp=Psi0_Em(:,ikr)
					Ex0_Em(ikr)=Ex0_Em(ikpr)
					Psi0_Em(:,ikr)=Psi0_Em(:,ikpr)
					Ex0_Em(ikpr)=tmpr
					Psi0_Em(:,ikpr)=Psi_tmp
				end if 
				
				if (Ex1_Em(ikr) .lt. Ex1_Em(ikpr)) then
					tmpr=Ex1_Em(ikr)
					Psi_tmp=Psi1_Em(:,ikr)
					Ex1_Em(ikr)=Ex1_Em(ikpr)
					Psi1_Em(:,ikr)=Psi1_Em(:,ikpr)
					Ex1_Em(ikpr)=tmpr
					Psi1_Em(:,ikpr)=Psi_tmp
				end if 	
			end do
		end do
		return
	end subroutine calculate_exciton_energy	

end module E_exciton_energy_mod