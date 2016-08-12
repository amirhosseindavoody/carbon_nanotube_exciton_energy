module cnt_geometry_mod
	implicit none
	private
	public :: cnt_geometry

contains
	!**********************************************************************************************************************
	! This subroutines interprets the input arguments of the simulation
	!**********************************************************************************************************************

	subroutine cnt_geometry()
		use comparams, only: currcnt
		use math_functions_mod, only: gcd, my_norm2
		use physical_constant_mod, only: a_l, pi
		use write_log_mod, only: writeLog

		integer :: tmpi,i,j,k
		integer :: MC
		integer :: dR = 1
        integer :: t1,t2
		real*8 :: p1,p2,p,q
		real*8 :: cosTh, sinTh
		real*8, dimension(2,2) :: Rot
		character(len=1000) :: logInput

		! unit vectors and reciprocal lattice vectors************************************************************************
		currcnt%a1=(/dsqrt(3.d0)/2.d0*a_l, +1.d0/2.d0*a_l/)
		currcnt%a2=(/dsqrt(3.d0)/2.d0*a_l, -1.d0/2.d0*a_l/)
		currcnt%b1=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, +1.d0*2.d0*pi/a_l/)
		currcnt%b2=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, -1.d0*2.d0*pi/a_l/)

		currcnt%aCC_vec=1.d0/3.d0*(currcnt%a1+currcnt%a2)

		! calculate chirality and translational vectors of CNT unit cell.****************************************************
		currcnt%ch_vec=dble(currcnt%n_ch)*currcnt%a1+dble(currcnt%m_ch)*currcnt%a2

		currcnt%len_ch=a_l*dsqrt(dble(currcnt%n_ch)**2+dble(currcnt%m_ch)**2+dble(currcnt%n_ch)*dble(currcnt%m_ch))
		currcnt%radius=currcnt%len_ch/2.d0/pi

		call gcd(dR,2*currcnt%n_ch+currcnt%m_ch,2*currcnt%m_ch+currcnt%n_ch)

		t1=int((2.d0*dble(currcnt%m_ch)+dble(currcnt%n_ch))/dble(dR))
		t2=-int((2.d0*dble(currcnt%n_ch)+dble(currcnt%m_ch))/dble(dR))

		currcnt%t_vec=dble(t1)*currcnt%a1+ dble(t2)*currcnt%a2

		currcnt%Nu=2*(currcnt%n_ch**2+currcnt%m_ch**2+currcnt%n_ch*currcnt%m_ch)/dR

		p1=max((dble(currcnt%Nu)/dble(currcnt%n_ch)+1.d0/dble(t1))/(dble(currcnt%m_ch)/dble(currcnt%n_ch)- dble(t2)/dble(t1)),(1.d0/dble(currcnt%n_ch)+1.d0/dble(t1))/(dble(currcnt%m_ch)/dble(currcnt%n_ch)-dble(t2)/dble(t1)))
		p2=min((dble(currcnt%Nu)/dble(currcnt%n_ch)+1.d0/dble(t1))/(dble(currcnt%m_ch)/dble(currcnt%n_ch)- dble(t2)/dble(t1)),(1.d0/dble(currcnt%n_ch)+1.d0/dble(t1))/(dble(currcnt%m_ch)/dble(currcnt%n_ch)-dble(t2)/dble(t1)))

		do i=ceiling(p2),floor(p1)
			p=dble(i)
			q=dble(t2)/dble(t1)*p+1.d0/dble(t1)
			if (q .eq. dble(ceiling(q))) then
				exit
			else if (i==floor(p1)) then
				call writeLog("MC not found")
				call exit()
			end if
		end do

		call gcd(tmpi,int(abs(p)),int(abs(q)))
		p=p/dble(tmpi)
		q=q/dble(tmpi)
		MC=int(dble(currcnt%m_ch)*p-dble(currcnt%n_ch)*q)

		! rotate basis vectors so that ch_vec is along x-axis
		cosTh=currcnt%ch_vec(1)/my_norm2(currcnt%ch_vec)
		sinTh=currcnt%ch_vec(2)/my_norm2(currcnt%ch_vec)
		Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
		currcnt%ch_vec=matmul(Rot,currcnt%ch_vec)
		currcnt%t_vec=matmul(Rot,currcnt%t_vec)
		currcnt%a1=matmul(Rot,currcnt%a1)
		currcnt%a2=matmul(Rot,currcnt%a2)
		currcnt%b1=matmul(Rot,currcnt%b1)
		currcnt%b2=matmul(Rot,currcnt%b2)
		currcnt%aCC_vec=matmul(Rot,currcnt%aCC_vec)

		! calculate reciprocal lattice of CNT.*******************************************************************************
		currcnt%dk=my_norm2(currcnt%b1)/(dble(currcnt%nkg)-1.d0)
		currcnt%dkx = currcnt%dk/currcnt%dk_dkx_ratio
		currcnt%K1=(- t2*currcnt%b1+ dble(t1)*currcnt%b2)/(dble(currcnt%Nu))
		currcnt%K2=(dble(currcnt%m_ch)*currcnt%b1-dble(currcnt%n_ch)*currcnt%b2)/dble(currcnt%Nu)
		currcnt%K2=currcnt%K2/my_norm2(currcnt%K2)

		! calculate coordinates of atoms in the unwarped CNT unit cell.******************************************************
		allocate(currcnt%posA(currcnt%Nu,2))
		allocate(currcnt%posB(currcnt%Nu,2))

		open(unit=100,file='posA.dat',status="unknown")
		open(unit=101,file='posB.dat',status="unknown")

		k=0
		do i=0, t1+currcnt%n_ch
			do j= t2,currcnt%m_ch
				if ((dble(t2*i)/dble(t1) .le. dble(j)) .and. (dble(currcnt%m_ch*i)/dble(currcnt%n_ch) .ge. dble(j)) .and. (dble( t2)/dble(t1)*dble(i-currcnt%n_ch) .gt. dble(j-currcnt%m_ch)) .and. (dble(currcnt%m_ch)/dble(currcnt%n_ch)*dble(i-t1) .lt. dble(j-t2))) then
					k=k+1
					currcnt%posA(k,1)=dble(i)*currcnt%a1(1)+dble(j)*currcnt%a2(1)
					currcnt%posA(k,2)=dble(i)*currcnt%a1(2)+dble(j)*currcnt%a2(2)
					currcnt%posB(k,1)=currcnt%posA(k,1)+currcnt%aCC_vec(1)
					currcnt%posB(k,2)=currcnt%posA(k,2)+currcnt%aCC_vec(2)

					if (currcnt%posA(k,1) .gt. currcnt%ch_vec(1)) currcnt%posA(k,1)=currcnt%posA(k,1)-currcnt%ch_vec(1);
					if (currcnt%posA(k,1) .lt. 0) currcnt%posA(k,1)=currcnt%posA(k,1)+currcnt%ch_vec(1);
					if (currcnt%posA(k,2) .gt. currcnt%t_vec(2)) currcnt%posA(k,2)=currcnt%posA(k,2)-currcnt%t_vec(2);
					if (currcnt%posA(k,2) .lt. 0) currcnt%posA(k,2)=currcnt%posA(k,2)+currcnt%t_vec(2);

					if (currcnt%posB(k,1) .gt. currcnt%ch_vec(1)) currcnt%posB(k,1)=currcnt%posB(k,1)-currcnt%ch_vec(1);
					if (currcnt%posB(k,1) .lt. 0) currcnt%posB(k,1)=currcnt%posB(k,1)+currcnt%ch_vec(1);
					if (currcnt%posB(k,2) .gt. currcnt%t_vec(2)) currcnt%posB(k,2)=currcnt%posB(k,2)-currcnt%t_vec(2);
					if (currcnt%posB(k,2) .lt. 0) currcnt%posB(k,2)=currcnt%posB(k,2)+currcnt%t_vec(2);

					write(100,'(E16.8,E16.8)') currcnt%posA(k,1), currcnt%posA(k,2)
					write(101,'(E16.8,E16.8)') currcnt%posB(k,1), currcnt%posB(k,2)

				endif
			enddo
		enddo

		close(100)
		close(101)

		if (k .ne. currcnt%Nu) then
			call writeLog("*** Error in calculating atom positions ***")
			call exit()
		end if

		! calculate distances between atoms in a warped CNT unit cell.*******************************************************
		allocate(currcnt%posAA(currcnt%Nu,2))
		allocate(currcnt%posAB(currcnt%Nu,2))
		allocate(currcnt%posBA(currcnt%Nu,2))
		allocate(currcnt%posBB(currcnt%Nu,2))

		do i=1,currcnt%Nu
			currcnt%posAA(i,:)=currcnt%posA(i,:)-currcnt%posA(1,:)
			currcnt%posAB(i,:)=currcnt%posA(i,:)-currcnt%posB(1,:)
			currcnt%posBA(i,:)=currcnt%posB(i,:)-currcnt%posA(1,:)
			currcnt%posBB(i,:)=currcnt%posB(i,:)-currcnt%posB(1,:)
			if (currcnt%posAA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAA(i,1)=currcnt%posAA(i,1)-currcnt%ch_vec(1)
			if (currcnt%posAB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posAB(i,1)=currcnt%posAB(i,1)-currcnt%ch_vec(1)
			if (currcnt%posBA(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBA(i,1)=currcnt%posBA(i,1)-currcnt%ch_vec(1)
			if (currcnt%posBB(i,1) .gt. currcnt%ch_vec(1)/2.d0) currcnt%posBB(i,1)=currcnt%posBB(i,1)-currcnt%ch_vec(1)

! 			write(100,'(E16.8,E16.8)') currcnt%posAA(k,1), currcnt%posAA(k,2)
! 			write(101,'(E16.8,E16.8)') currcnt%posBA(k,1), currcnt%posBA(k,2)

		end do

		! write down important informations into the output file.************************************************************
		write(logInput,'(A, I0, A, I0, A)') "chirality = (",currcnt%n_ch, " , ", currcnt%m_ch, ")"
		call writeLog(trim(logInput))
		write(logInput,'(A, E16.8)') "radius = ", currcnt%radius
		call writeLog(trim(logInput))

		call writeLog(new_line('A')//"geometrical properties **********************************")
		write(logInput,'(SP, A, E9.2, E9.2)') "a1=",currcnt%a1(1), currcnt%a1(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "a2=",currcnt%a2(1), currcnt%a2(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "b1=",currcnt%b1(1), currcnt%b1(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "b2=",currcnt%b2(1), currcnt%b2(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "aCC_vec=",currcnt%aCC_vec(1), currcnt%aCC_vec(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "ch_vec=",currcnt%ch_vec(1), currcnt%ch_vec(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2, E9.2)') "t_vec=",currcnt%t_vec(1), currcnt%t_vec(2)
		call writeLog(trim(logInput))
		write(logInput,'(SP, A, E9.2)') "len_ch=",currcnt%len_ch
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "Nu=",currcnt%Nu
		call writeLog(trim(logInput))
		write(logInput,'(A, I0)') "MC=",MC
		call writeLog(trim(logInput))
		write(logInput,'(A, E9.2)') "dk= ",currcnt%dk
		call writeLog(trim(logInput))
		write(logInput,'(A, E9.2)') "dkx=",currcnt%dkx
		call writeLog(trim(logInput))

		return
	end subroutine cnt_geometry

end module cnt_geometry_mod
