module dielectric_fn_mod
	implicit none
	private
	public :: calculate_dielectric_fn

contains
	subroutine calculate_dielectric_fn()
		use cnt_band_structure_mod, only: grapheneEnergy
		use comparams, only: currcnt, flg_dielectric
		use math_functions_mod, only: my_norm2
		use physical_constant_mod, only: i1, pi, q0, Upp, eps0
		use write_log_mod, only: writeLog
		
		integer :: ik, mu_k, iq, mu_q, ikq, mu_kq
		integer :: u,b
		real*8 :: tmpr
		real*8, dimension(2) :: deltaR
		real*8, dimension(2) :: kq
		real*8, dimension(2) :: E_tmp
		real*8, dimension(:,:), allocatable :: PI_q,v_q
		real*8, dimension(:,:,:), allocatable :: E_k
		real*8, dimension(:,:,:), allocatable :: q_vec
		complex*16, dimension(:,:,:), allocatable :: Cc_k, Cv_k
		complex*16, dimension(2) :: Cc_tmp, Cv_tmp
		character(len=200) :: logInput
		
		! initialize variables.**********************************************************************************************
		allocate(PI_q(1-currcnt%Nu:currcnt%Nu-1,currcnt%iq_min:currcnt%iq_max))
		allocate(v_q(1-currcnt%Nu:currcnt%Nu-1,currcnt%iq_min:currcnt%iq_max))
		allocate(q_vec(1-currcnt%Nu:currcnt%Nu-1,currcnt%iq_min:currcnt%iq_max,2))
		allocate(currcnt%eps_q(1-currcnt%Nu:currcnt%Nu-1,currcnt%iq_min:currcnt%iq_max))
		allocate(currcnt%v_FT(1-currcnt%Nu:currcnt%Nu-1,currcnt%iq_min:currcnt%iq_max,2,2))
		
		allocate(Cc_k(2-3*currcnt%Nu/2:3*currcnt%Nu/2-1,currcnt%iq_min+currcnt%ikc_min:currcnt%iq_max+currcnt%ikc_max,2))
		allocate(Cv_k(2-3*currcnt%Nu/2:3*currcnt%Nu/2-1,currcnt%iq_min+currcnt%ikc_min:currcnt%iq_max+currcnt%ikc_max,2))
		allocate(E_k(2-3*currcnt%Nu/2:3*currcnt%Nu/2-1,currcnt%iq_min+currcnt%ikc_min:currcnt%iq_max+currcnt%ikc_max,2))
		
		do mu_q=1-currcnt%Nu,currcnt%Nu-1
			do iq=currcnt%iq_min,currcnt%iq_max
				PI_q(mu_q,iq)=0.d0
				v_q(mu_q,iq)=0.d0
				currcnt%eps_q(mu_q,iq)=0.d0
				currcnt%v_FT(mu_q,iq,1,1)=(0.d0,0.d0)
				currcnt%v_FT(mu_q,iq,1,2)=(0.d0,0.d0)
				currcnt%v_FT(mu_q,iq,2,1)=(0.d0,0.d0)
				currcnt%v_FT(mu_q,iq,2,2)=(0.d0,0.d0)
			enddo
		enddo
		
		if (flg_dielectric .eqv. .true.) then
		! calculate PI(q).***************************************************************************************************
			do mu_kq=2-3*currcnt%Nu/2,3*currcnt%Nu/2-1
				do ikq=(currcnt%iq_min+currcnt%ikc_min),(currcnt%iq_max+currcnt%ikc_max)
					kq=dble(mu_kq)*currcnt%K1+dble(ikq)*currcnt%dk*currcnt%K2
					call grapheneEnergy(currcnt,E_tmp,Cc_tmp,Cv_tmp,kq)
					E_k(mu_kq,ikq,:)=E_tmp
					Cc_k(mu_kq,ikq,:)=Cc_tmp
					Cv_k(mu_kq,ikq,:)=Cv_tmp
				enddo
			enddo
			
			call writeLog(new_line('A')//"Calculating PI_q ************************")

			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				write(logInput,*) "mu_q=", mu_q
				call writeLog(trim(logInput))
				do iq=currcnt%iq_min,currcnt%iq_max
					do mu_k=1-currcnt%Nu/2,currcnt%Nu/2
						do ik=currcnt%ikc_min,currcnt%ikc_max
							ikq=ik+iq
							mu_kq=mu_k+mu_q
							PI_q(mu_q,iq)=PI_q(mu_q,iq)+ &
								(abs (dot_product(Cv_k(mu_k,ik,:),Cc_k(mu_kq,ikq,:))))**2/(E_k(mu_kq,ikq,1)-E_k(mu_k,ik,2))+ &
								(abs (dot_product(Cc_k(mu_k,ik,:),Cv_k(mu_kq,ikq,:))))**2/(E_k(mu_k,ik,1)-E_k(mu_kq,ikq,2))
						enddo
					enddo
				enddo
			enddo
		
		! calculate v_FT and eps_q.******************************************************************************************
			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				do iq=currcnt%iq_min,currcnt%iq_max
					q_vec(mu_q,iq,1)=dble(mu_q)*currcnt%K1(1)+dble(iq)*currcnt%dk*currcnt%K2(1)
					q_vec(mu_q,iq,2)=dble(mu_q)*currcnt%K1(2)+dble(iq)*currcnt%dk*currcnt%K2(2)
				enddo
			enddo
			
			call writeLog(new_line('A')//"Calculating v_FT ************************")

			tmpr = 4.d0*pi*eps0/(q0**2)*Upp
			do u=-currcnt%nr,currcnt%nr
				write(logInput,*) "u=",u
				call writeLog(trim(logInput))
				do b=1,currcnt%Nu 
					deltaR=dble(u)*currcnt%t_vec+currcnt%posAA(b,:)
					currcnt%v_FT(:,:,1,1)=currcnt%v_FT(:,:,1,1)+exp(i1*dcmplx(q_vec(:,:,1)*deltaR(1)+q_vec(:,:,2)*deltaR(2)))* &
							dcmplx(Upp/sqrt((tmpr*my_norm2(deltaR))**2+1.d0))
					
					deltaR=dble(u)*currcnt%t_vec+currcnt%posBA(b,:)
					currcnt%v_FT(:,:,1,2)=currcnt%v_FT(:,:,1,2)+exp(i1*dcmplx(q_vec(:,:,1)*deltaR(1)+q_vec(:,:,2)*deltaR(2)))* &
							dcmplx(Upp/sqrt((tmpr*my_norm2(deltaR))**2+1.d0))
				enddo
			enddo
			
			currcnt%v_FT(:,:,2,2)=currcnt%v_FT(:,:,1,1)
			
			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				do iq=currcnt%iq_min,currcnt%iq_max
					currcnt%v_FT(mu_q,iq,2,1)=currcnt%v_FT(-mu_q,-iq,1,2)
				enddo
			enddo
			
			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				do iq=currcnt%iq_min,currcnt%iq_max
					v_q(mu_q,iq)=(1.d0/4.d0)*dreal(sum(currcnt%v_FT(mu_q,iq,:,:)))
				enddo
			enddo
			

			currcnt%v_FT=currcnt%v_FT/dcmplx((currcnt%ikc_max-currcnt%ikc_min+1)*currcnt%Nu) ! scaling of v_FT is done to make interaction kernel elements smaller so as the system grows in size the eigen values do not scale. (then number of k points is a function of the length of cnt)
			v_q=v_q/dble((currcnt%ikc_max-currcnt%ikc_min+1)*currcnt%Nu) ! the scaling of v_q is done to make eps_q almost size independent.
			
			currcnt%eps_q=1.d0+v_q*PI_q
			
			
		! save eps_q, and v_FT.**********************************************************************************************
			
			open(unit=100,file='qVec.dat',status="unknown")
			do iq=currcnt%iq_min,currcnt%iq_max
				write(100,'(E16.8)', advance='no') dble(iq)*currcnt%dk
			end do
			close(100)

			open(unit=100,file='PI_q.dat',status="unknown")
			open(unit=101,file='v_FT.dat',status="unknown")
			open(unit=102,file='eps_q.dat',status="unknown")

			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				do iq=currcnt%iq_min,currcnt%iq_max
					write(100,'(E16.8)', advance='no') PI_q(mu_q,iq)
					write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,1,1) 
					write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,1,2) 
					write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,2,1) 
					write(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,2,2)
					write(102,'(E16.8)', advance='no') currcnt%eps_q(mu_q,iq) 
				end do
				
				write(100,*)
				write(101,*)
				write(102,*)
			enddo

			close(100)
			close(101)
			close(102)
		
		else

			open(unit=100,file='PI_q.dat',status="old")
			open(unit=101,file='v_FT.dat',status="old")
			open(unit=102,file='eps_q.dat',status="old")

			do mu_q=1-currcnt%Nu,currcnt%Nu-1
				do iq=currcnt%iq_min,currcnt%iq_max
					read(100,'(E16.8)', advance='no') PI_q(mu_q,iq)
					read(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,1,1) 
					read(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,1,2) 
					read(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,2,1) 
					read(101,'(E16.8,E16.8)', advance='no') currcnt%v_FT(mu_q,iq,2,2)
					read(102,'(E16.8)', advance='no') currcnt%eps_q(mu_q,iq) 
				end do
				
				read(100,*)
				read(101,*)
				read(102,*)
			enddo

			close(100)
			close(101)
			close(102)

		endif

	end subroutine calculate_dielectric_fn

end module dielectric_fn_mod