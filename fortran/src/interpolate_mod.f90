module interpolate_mod
	implicit none
	private
	public :: interpolate_energy, interpolate_dielectric_fn

contains

	!**********************************************************************************************************************
	! This subroutines interpolates the calculated self-energy and tight-binding energies to get a finer mesh results
	!**********************************************************************************************************************

	subroutine interpolate_energy()
		use cnt_band_structure_mod, only: grapheneEnergy
		use comparams, only: currcnt
		use math_functions_mod, only: polint
		use write_log_mod, only: writeLog
		
		integer :: ik_low_fine, ik_high_fine
		integer :: ik, mu_k
		integer :: tmpi, i, j, ik_tmp
		integer, parameter :: nInterpolate=11
		real*8 :: error, x
		real*8, dimension(2) :: k
		real*8, dimension(2) :: E_tmp
		real*8, dimension(nInterpolate) :: ya, xa
		complex*16, dimension(2) :: Cc_tmp,Cv_tmp
		character(len=200) :: logInput
		
		ik_low_fine = currcnt%ik_low * currcnt%dk_dkx_ratio
		ik_high_fine = currcnt%ik_high * currcnt%dk_dkx_ratio

		! initialize variables.**********************************************************************************************
		allocate(currcnt%Ek_fine(2,ik_low_fine:ik_high_fine,2))
		allocate(currcnt%Sk_fine(2,ik_low_fine:ik_high_fine,2))
		allocate(currcnt%Cc_fine(2,ik_low_fine:ik_high_fine,2))
		allocate(currcnt%Cv_fine(2,ik_low_fine:ik_high_fine,2))
		
		call writeLog(new_line('A')//"Interpolating self energy **************************")

		! calculate the self energies and the tight-binding energies and coefficients.***************************************
		do ik=ik_low_fine,ik_high_fine
			
! 			write(logInput,'(A, I0)')'ik=',ik
! 			call writeLog(trim(logInput))

			if (mod(ik,currcnt%dk_dkx_ratio) .eq. 0) then
				currcnt%Ek_fine(:,ik,:) = currcnt%Ek(:,ik/(currcnt%dk_dkx_ratio),:)
				currcnt%Sk_fine(:,ik,:) = currcnt%Sk(:,ik/(currcnt%dk_dkx_ratio),:)
				currcnt%Cc_fine(:,ik,:) = currcnt%Cc(:,ik/(currcnt%dk_dkx_ratio),:)
				currcnt%Cv_fine(:,ik,:) = currcnt%Cv(:,ik/(currcnt%dk_dkx_ratio),:)
			else

				mu_k=currcnt%min_sub(currcnt%i_sub)
				k=dble(mu_k)*currcnt%K1+dble(ik)*currcnt%dkx*currcnt%K2
				call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
				currcnt%Ek_fine(1,ik,:)=E_tmp
				currcnt%Cc_fine(1,ik,:)=Cc_tmp
				currcnt%Cv_fine(1,ik,:)=Cv_tmp

				mu_k=-currcnt%min_sub(currcnt%i_sub)
				k=dble(mu_k)*currcnt%K1+dble(ik)*currcnt%dkx*currcnt%K2
				call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
				currcnt%Ek_fine(2,ik,:)=E_tmp
				currcnt%Cc_fine(2,ik,:)=Cc_tmp
				currcnt%Cv_fine(2,ik,:)=Cv_tmp				
				
				x = dble(ik)*currcnt%dkx

				tmpi = nint(dble(ik)/dble(currcnt%dk_dkx_ratio))
				if (tmpi .gt. (currcnt%ik_high-(nInterpolate-1)/2)) then
					do i=1,2
						do j=1,2
							do ik_tmp=1,nInterpolate
								xa(ik_tmp) = dble(currcnt%ik_high-nInterpolate+ik_tmp)*currcnt%dk
								ya(ik_tmp) = currcnt%Sk(i,currcnt%ik_high-nInterpolate+ik_tmp,j)
							enddo
							call polint(xa,ya,nInterpolate,x,currcnt%Sk_fine(i,ik,j),error)
						enddo
					enddo
				elseif ( tmpi .lt. (currcnt%ik_low+(nInterpolate-1)/2)) then
					do i=1,2
						do j=1,2
							do ik_tmp=1,nInterpolate
								xa(ik_tmp) = dble(currcnt%ik_low+ik_tmp-1)*currcnt%dk
								ya(ik_tmp) = currcnt%Sk(i,currcnt%ik_low+ik_tmp-1,j)
							enddo
							call polint(xa,ya,nInterpolate,x,currcnt%Sk_fine(i,ik,j),error)
						enddo
					enddo
				else
					do i=1,2
						do j=1,2
							do ik_tmp=-(nInterpolate-1)/2,(nInterpolate-1)/2
								xa(ik_tmp+(nInterpolate+1)/2) = dble(tmpi+ik_tmp)*currcnt%dk
								ya(ik_tmp+(nInterpolate+1)/2) = currcnt%Sk(i,tmpi+ik_tmp,j)
							enddo
							
							call polint(xa,ya,nInterpolate,x,currcnt%Sk_fine(i,ik,j),error)
						enddo
					enddo
				endif

			end if
		enddo
		
		!save the results.***************************************************************************************************
		open(unit=100,file='CondSelfEnergy_fine.dat',status="unknown")
		open(unit=101,file='ValeSelfEnergy_fine.dat',status="unknown")
		open(unit=102,file='CondBand_fine.dat',status="unknown")
		open(unit=103,file='ValeBand_fine.dat',status="unknown")
		open(unit=104,file='kVec_fine.dat',status="unknown")

		do ik=ik_low_fine,ik_high_fine
			write(100,'(E16.8)', advance='no') currcnt%Sk_fine(1,ik,1) 
			write(101,'(E16.8)', advance='no') currcnt%Sk_fine(1,ik,2)
			write(102,'(E16.8)', advance='no') currcnt%Ek_fine(1,ik,1) 
			write(103,'(E16.8)', advance='no') currcnt%Ek_fine(1,ik,2)
			write(104,'(E16.8)', advance='no') dble(ik)*currcnt%dkx
		enddo
		
		write(100,*)
		write(101,*)
		write(102,*)
		write(103,*)
		
		do ik=ik_low_fine,ik_high_fine
			write(100,'(E16.8)', advance='no') currcnt%Sk_fine(2,ik,1) 
			write(101,'(E16.8)', advance='no') currcnt%Sk_fine(2,ik,2) 
			write(102,'(E16.8)', advance='no') currcnt%Ek_fine(2,ik,1) 
			write(103,'(E16.8)', advance='no') currcnt%Ek_fine(2,ik,2) 
		enddo

		close(100)
		close(101)
		close(102)
		close(103)
		close(104)
		
	end subroutine interpolate_energy

	!**********************************************************************************************************************
	! This subroutines interpolates the calculated v_FT and eps_q to a finer mesh size
	!**********************************************************************************************************************

	subroutine interpolate_dielectric_fn()
		use comparams, only: currcnt
		use math_functions_mod, only: polint
		use write_log_mod, only: writeLog
		
		integer :: iq_min_fine, iq_max_fine
		integer :: iq, mu_q
		integer :: iq_tmp, tmpi
		integer :: i, j
		integer, parameter :: nInterpolate=11
		real*8 :: x, error, tmpr1, tmpr2
		real*8, dimension(nInterpolate) :: ya, xa, ya1, ya2
		character(len=200) :: logInput
		
		iq_min_fine = currcnt%iq_min * currcnt%dk_dkx_ratio
		iq_max_fine = currcnt%iq_max * currcnt%dk_dkx_ratio

		! initialize variables.**********************************************************************************************
		allocate(currcnt%eps_q_fine(1-currcnt%Nu:currcnt%Nu-1,iq_min_fine:iq_max_fine))
		allocate(currcnt%v_FT_fine(1-currcnt%Nu:currcnt%Nu-1,iq_min_fine:iq_max_fine,2,2))
		
		call writeLog(new_line('A')//"Interpolating dielectric function **************************")

		! Interpolate eps_q and v_FT
		do iq=iq_min_fine,iq_max_fine
			
! 			write(logInput,'(A, I0)')'iq=',iq
! 			call writeLog(trim(logInput))

			if (mod(iq,currcnt%dk_dkx_ratio) .eq. 0) then
				currcnt%eps_q_fine(:,iq) = currcnt%eps_q(:,iq/(currcnt%dk_dkx_ratio))
				currcnt%v_FT_fine(:,iq,:,:) = currcnt%v_FT(:,iq/(currcnt%dk_dkx_ratio),:,:)
			else
				do mu_q=1-currcnt%Nu,currcnt%Nu-1
					x = dble(iq)*currcnt%dkx
					tmpi = nint(dble(iq)/dble(currcnt%dk_dkx_ratio))

					if (tmpi .gt. (currcnt%iq_max-(nInterpolate-1)/2)) then
						
						do iq_tmp=1,nInterpolate
							xa(iq_tmp) = dble(currcnt%iq_max-nInterpolate+iq_tmp)*currcnt%dk
							ya(iq_tmp) = currcnt%eps_q(mu_q,currcnt%iq_max-nInterpolate+iq_tmp)
						enddo
						call polint(xa,ya,nInterpolate,x,currcnt%eps_q_fine(mu_q,iq),error)

						do i=1,2
							do j=1,2
								do iq_tmp=1,nInterpolate
									xa(iq_tmp) = dble(currcnt%iq_max-nInterpolate+iq_tmp)*currcnt%dk
									ya1(iq_tmp) = real(currcnt%v_FT(mu_q,currcnt%iq_max-nInterpolate+iq_tmp,i,j))
									ya2(iq_tmp) = aimag(currcnt%v_FT(mu_q,currcnt%iq_max-nInterpolate+iq_tmp,i,j))
								enddo
								call polint(xa,ya1,nInterpolate,x,tmpr1,error)
								call polint(xa,ya2,nInterpolate,x,tmpr2,error)
								currcnt%v_FT_fine(mu_q,iq,i,j) = complex(tmpr1,tmpr2)
							enddo
						enddo

					elseif ( tmpi .lt. (currcnt%iq_min+(nInterpolate-1)/2)) then

						do iq_tmp=1,nInterpolate
							xa(iq_tmp) = dble(currcnt%iq_min+iq_tmp-1)*currcnt%dk
							ya(iq_tmp) = currcnt%eps_q(mu_q,currcnt%iq_min+iq_tmp-1)
						enddo
						call polint(xa,ya,nInterpolate,x,currcnt%eps_q_fine(mu_q,iq),error)

						do i=1,2
							do j=1,2
								do iq_tmp=1,nInterpolate
									xa(iq_tmp) = dble(currcnt%iq_min+iq_tmp-1)*currcnt%dk
									ya1(iq_tmp) = real(currcnt%v_FT(mu_q,currcnt%iq_min+iq_tmp-1,i,j))
									ya2(iq_tmp) = aimag(currcnt%v_FT(mu_q,currcnt%iq_min+iq_tmp-1,i,j))
								enddo
								call polint(xa,ya1,nInterpolate,x,tmpr1,error)
								call polint(xa,ya2,nInterpolate,x,tmpr2,error)
								currcnt%v_FT_fine(mu_q,iq,i,j) = complex(tmpr1,tmpr2)
							enddo
						enddo


					else
						do iq_tmp=-(nInterpolate-1)/2,(nInterpolate-1)/2
							xa(iq_tmp+(nInterpolate+1)/2) = dble(tmpi+iq_tmp)*currcnt%dk
							ya(iq_tmp+(nInterpolate+1)/2) = currcnt%eps_q(mu_q,tmpi+iq_tmp)
						enddo
						call polint(xa,ya,nInterpolate,x,currcnt%eps_q_fine(mu_q,iq),error)

						do i=1,2
							do j=1,2
								do iq_tmp=-(nInterpolate-1)/2,(nInterpolate-1)/2
									xa(iq_tmp+(nInterpolate+1)/2) = dble(tmpi+iq_tmp)*currcnt%dk
									ya1(iq_tmp+(nInterpolate+1)/2) = real(currcnt%v_FT(mu_q,tmpi+iq_tmp,i,j))
									ya2(iq_tmp+(nInterpolate+1)/2) = aimag(currcnt%v_FT(mu_q,tmpi+iq_tmp,i,j))
								enddo
								call polint(xa,ya1,nInterpolate,x,tmpr1,error)
								call polint(xa,ya2,nInterpolate,x,tmpr2,error)
								currcnt%v_FT_fine(mu_q,iq,i,j) = complex(tmpr1,tmpr2)
							enddo
						enddo


					endif
				enddo			
				
				

				

			end if
		enddo
			
			
		! save eps_q, and v_FT.**********************************************************************************************
			
		open(unit=100,file='qVec_fine.dat',status="unknown")
		do iq=iq_min_fine,iq_max_fine
			write(100,'(E16.8)', advance='no') dble(iq)*currcnt%dkx
		enddo
		close(100)

		open(unit=101,file='v_FT_fine.dat',status="unknown")
		open(unit=102,file='eps_q_fine.dat',status="unknown")
	
		do mu_q=1-currcnt%Nu,currcnt%Nu-1
			do iq=iq_min_fine,iq_max_fine
				write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT_fine(mu_q,iq,1,1) 
				write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT_fine(mu_q,iq,1,2) 
				write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT_fine(mu_q,iq,2,1) 
				write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT_fine(mu_q,iq,2,2)
				write(102,'(E16.8)', advance='no') currcnt%eps_q_fine(mu_q,iq) 
			end do
				
			write(101,*)
			write(102,*)
		enddo

		close(101)
		close(102)

	end subroutine interpolate_dielectric_fn

end module interpolate_mod