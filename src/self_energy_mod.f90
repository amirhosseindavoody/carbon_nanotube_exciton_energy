module self_energy_mod
	implicit none
	private
	public :: calculate_self_energy

contains	
	subroutine calculate_self_energy()
		use cnt_band_structure_mod, only: grapheneEnergy
		use comparams, only: currcnt
		use write_log_mod, only: writeLog
		
		integer :: ik, mu_k, iq, mu_q
		real*8, dimension(2) :: E_tmp
		real*8, dimension(2) :: k,kq
		real*8, dimension(2) :: E_kq
		complex*16, dimension(2) :: Cc_kq,Cv_kq
		complex*16, dimension(2) :: Cc_tmp,Cv_tmp
		character(len=200) :: logInput
		
		! initialize variables.**********************************************************************************************
		allocate(currcnt%Ek(2,currcnt%ik_low:currcnt%ik_high,2))
		allocate(currcnt%Sk(2,currcnt%ik_low:currcnt%ik_high,2))
		allocate(currcnt%Cc(2,currcnt%ik_low:currcnt%ik_high,2))
		allocate(currcnt%Cv(2,currcnt%ik_low:currcnt%ik_high,2))
		
		do ik=currcnt%ik_low,currcnt%ik_high
			currcnt%Sk(1,ik,1)=0.d0
			currcnt%Sk(1,ik,2)=0.d0
			currcnt%Sk(2,ik,1)=0.d0
			currcnt%Sk(2,ik,2)=0.d0
		enddo
		
		call writeLog(new_line('A')//"Calculating self energy **************************")

		! calculate the self energies and the tight-binding energies and coefficients.***************************************
		do ik=currcnt%ik_low,currcnt%ik_high
			
			write(logInput,*)'ik=',ik
			call writeLog(trim(logInput))
			
			mu_k=currcnt%min_sub(currcnt%i_sub)
			k=dble(mu_k)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
			currcnt%Ek(1,ik,:)=E_tmp
			currcnt%Cc(1,ik,:)=Cc_tmp
			currcnt%Cv(1,ik,:)=Cv_tmp
			
			!calculate self energy of the first subband
			do mu_q=1-currcnt%Nu/2,currcnt%Nu/2
				do iq=currcnt%ikc_min,currcnt%ikc_max
					kq=dble(mu_k+mu_q)*currcnt%K1+dble(ik+iq)*currcnt%dk*currcnt%K2
					call grapheneEnergy(currcnt,E_kq,Cc_kq,Cv_kq,kq) 
					currcnt%Sk(1,ik,1)=currcnt%Sk(1,ik,1)-dreal(conjg(currcnt%Cc(1,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,1)*currcnt%Cc(1,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cc(1,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,2)*currcnt%Cc(1,ik,2)*conjg(Cv_kq(2))+ &
																conjg(currcnt%Cc(1,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,1)*currcnt%Cc(1,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cc(1,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,2)*currcnt%Cc(1,ik,2)*conjg(Cv_kq(2)))/(currcnt%kappa*currcnt%eps_q(-mu_q,-iq))
					
					currcnt%Sk(1,ik,2)=currcnt%Sk(1,ik,2)-dreal(conjg(currcnt%Cv(1,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,1)*currcnt%Cv(1,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cv(1,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,2)*currcnt%Cv(1,ik,2)*conjg(Cv_kq(2))+ &
																conjg(currcnt%Cv(1,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,1)*currcnt%Cv(1,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cv(1,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,2)*currcnt%Cv(1,ik,2)*conjg(Cv_kq(2)))/(currcnt%kappa*currcnt%eps_q(-mu_q,-iq))
				enddo
			enddo
			
			!calculate the self energy of the second subband
			mu_k=-currcnt%min_sub(currcnt%i_sub)
			k=dble(mu_k)*currcnt%K1+dble(ik)*currcnt%dk*currcnt%K2
			call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,k)
			currcnt%Ek(2,ik,:) = E_tmp
			currcnt%Cc(2,ik,:) = Cc_tmp
			currcnt%Cv(2,ik,:) = Cv_tmp
			
			do mu_q=1-currcnt%Nu/2,currcnt%Nu/2
				do iq=currcnt%ikc_min,currcnt%ikc_max
					kq=dble(mu_k+mu_q)*currcnt%K1+dble(ik+iq)*currcnt%dk*currcnt%K2
					call grapheneEnergy(currcnt,E_kq,Cc_kq,Cv_kq,kq) 
					currcnt%Sk(2,ik,1)=currcnt%Sk(2,ik,1)-dreal(conjg(currcnt%Cc(2,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,1)*currcnt%Cc(2,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cc(2,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,2)*currcnt%Cc(2,ik,2)*conjg(Cv_kq(2))+ &
																conjg(currcnt%Cc(2,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,1)*currcnt%Cc(2,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cc(2,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,2)*currcnt%Cc(2,ik,2)*conjg(Cv_kq(2)))/(currcnt%kappa*currcnt%eps_q(-mu_q,-iq))
					
					currcnt%Sk(2,ik,2)=currcnt%Sk(2,ik,2)-dreal(conjg(currcnt%Cv(2,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,1)*currcnt%Cv(2,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cv(2,ik,1))*Cv_kq(1)*currcnt%v_FT(-mu_q,-iq,1,2)*currcnt%Cv(2,ik,2)*conjg(Cv_kq(2))+ &
																conjg(currcnt%Cv(2,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,1)*currcnt%Cv(2,ik,1)*conjg(Cv_kq(1))+ &
																conjg(currcnt%Cv(2,ik,2))*Cv_kq(2)*currcnt%v_FT(-mu_q,-iq,2,2)*currcnt%Cv(2,ik,2)*conjg(Cv_kq(2)))/(currcnt%kappa*currcnt%eps_q(-mu_q,-iq))
				enddo
			enddo
		enddo
		
		!save the results.***************************************************************************************************
		open(unit=100,file='CondSelfEnergy_Sub.dat',status="unknown")
		open(unit=101,file='ValeSelfEnergy_Sub.dat',status="unknown")

		do ik=currcnt%ik_low,currcnt%ik_high
			write(100,'(E16.8)', advance='no') currcnt%Sk(1,ik,1) 
			write(101,'(E16.8)', advance='no') currcnt%Sk(1,ik,2) 
		enddo
		
		write(100,*)
		write(101,*)
		
		do ik=currcnt%ik_low,currcnt%ik_high
			write(100,'(E16.8)', advance='no') currcnt%Sk(2,ik,1) 
			write(101,'(E16.8)', advance='no') currcnt%Sk(2,ik,2) 
		enddo

		close(100)
		close(101)
		
	end subroutine calculate_self_energy

end module self_energy_mod