module math_functions_mod
	implicit none
	private
	public :: gcd, bessk0, eig, polint, my_norm2

contains
	!**********************************************************************************************************************
	! This subroutine calculates the greatest common divisor of the arguments na and nb
	!**********************************************************************************************************************
	
	subroutine gcd(ngcd,na,nb)
	integer, intent(in) :: na,nb
	integer, intent(out) :: ngcd
	integer :: ia,ib,itemp

	ia=na
	ib=nb
	do while (ib .ne. 0)
		itemp=ia
		ia=ib
		ib=mod(itemp,ib)
	end do
	ngcd=ia
	return
	end subroutine gcd
		
	!**********************************************************************************************************************
	! This function calculates the modified bessel function of the first kind with parameter nu=0: I0(x)
	!**********************************************************************************************************************
	real*8 function bessi0(x)
		real*8 :: x
		real*8 :: ax
		real*8 :: y
		real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7, q8, q9
			
		data p1, p2, p3, p4, p5, p6, p7 /1.d0, 3.5156229d0, 3.0899424d0, 1.2067492d0, 0.2659732d0, 0.360768d-1, 0.45813d-2/
		data q1, q2, q3, q4, q5, q6, q7, q8, q9 /0.39894228d0, 0.1328592d-1, 0.225319d-2, -0.157565d-2, 0.916281d-2, -0.2057706d-1, 0.2635537d-1, -0.1647633d-1, 0.392377d-2/
		
		if (abs(x) .lt. 3.75d0) then
			y = (x/3.75d0)**2
			bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
		else
			ax = abs(x)
			y = 3.75d0/ax
			bessi0 = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
		end if
		return
	end function bessi0

	!**********************************************************************************************************************
	! This function calculates the modified bessel function of the second kind with parameter nu=0: K0(x)
	!**********************************************************************************************************************
	real*8 function bessk0(x)
		real*8 :: x
		real*8 :: y
		real*8, save :: p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, q5, q6, q7
		
		data p1, p2, p3, p4, p5, p6, p7 /-0.57721566d0, 0.42278420d0, 0.23069756d0, 0.3488590d-1, 0.262698d-2, 0.10750d-3, 0.74d-5/
		data q1, q2, q3, q4, q5, q6, q7 /1.25331414d0, -0.7832358d-1, 0.2189568d-1, -0.1062446d-1, 0.587872d-2, -0.251540d-2, 0.53208d-3/

		if (x .le. 2.d0) then
			y = x*x/4.d0
			bessk0=(-log(x/2.d0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
		else
			y=2.d0/x
			bessk0=(exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
		end if
		return 
	end function bessk0

	!**********************************************************************************************************************
	! This subroutine calculates the eigen values and eigen vectors of matrix given the size of it (nl) and gives back eigen vectors in A.
	!**********************************************************************************************************************

	subroutine eig(nl,matrix,A,W)

		interface 
			subroutine ZHEEV( JOBZ, UPLO,	N, A, LDA, W, WORK, LWORK, RWORK, INFO )
				character (len=1) :: JOBZ, UPLO
				integer :: INFO
				integer :: LDA, LWORK, N
				complex*16, dimension(N,N) :: A
				real*8, dimension(N) :: W
				complex*16, dimension(2*N-1) :: WORK
				real*8, dimension(3*N-2) :: RWORK
			end subroutine ZHEEV
		end interface

		integer, intent(in) :: nl
		complex*16, dimension(nl,nl), intent(inout) :: matrix
		complex*16, dimension(nl,nl), intent(out) :: A
		complex*16, dimension(2*nl-1) :: WORK
		character (len=1) :: JOBZ, UPLO
		integer :: INFO, LDA, LWORK, N
		real*8, dimension(nl) :: W
		real*8, dimension(3*nl-2) :: RWORK
		
		JOBZ = 'V'
		UPLO = 'L'
		N = size(matrix,1)
		A = matrix
		LDA = N
		LWORK = 2*size(matrix,1)-1
		
		call ZHEEV( JOBZ, UPLO,	N, A, LDA, W, WORK, LWORK, RWORK, INFO )
		
		if (INFO .ne. 0) then
			write(*,"(A50,I1.1)") "ERROR: calculation of eigen value failed , INFO = ", INFO
			stop
		end if
		
		return
	end subroutine eig

	!**********************************************************************************************************************
	! This subroutine interpolated arrays xa and ya using polynomial interpolation technique.
	! Given arrays xa and ya, each of length n, and given a value x, this routine returns a value y, and an error estimate.
	! if P(x) is athe polynomial of degree n-1 such that P(xa) = ya, for all values of arrays xa and ya, then the returned
	! value y = P(x).
	!**********************************************************************************************************************
	
	subroutine polint(xa, ya, n, x, y, dy)
		integer, intent(in) :: n
		real*8, intent(in) :: x
		real*8, intent(out) :: y, dy
		real*8, dimension(n), intent(in) :: xa, ya
		
		integer, parameter :: nmax=10
		integer :: i, m, ns
		real*8 :: den, dif, dift, ho, hp, w
		real*8, dimension(n) :: c, d

		ns = 1
		dif = abs(x-xa(1))

		do i=1,n
			dift = abs(x-xa(i))
			if (dift .lt. dif) then
				ns = i
				dif = dift
			endif
			c(i) = ya(i)
			d(i) = ya(i)
		enddo

		y = ya(ns)
		ns = ns-1

		do m=1,n-1
			do i=1, n-m
				ho=xa(i)-x
				hp=xa(i+m)-x
				w=c(i+1)-d(i)
				den=ho-hp
				if(den .eq. 0) then
					write(*,'(A)') "Faliure in polint"
					call exit()
				endif
				den = w/den
				d(i) = hp*den
				c(i) = ho*den
			enddo
			if (2*ns .lt. n-m) then
				dy = c(ns+1)
			else
				dy=d(ns)
				ns=ns-1
			endif
			y=y+dy
		enddo
		return
	end subroutine polint

	!**********************************************************************************************************************
	! This function calculates the magnitude of a 2D real*8 vector
	!**********************************************************************************************************************
	
	real*8 function my_norm2(my_vec)
		real*8, dimension(2), intent(in) :: my_vec
		
		my_norm2 = sqrt(my_vec(1)*my_vec(1)+my_vec(2)*my_vec(2))

		return 
	end function my_norm2

end module math_Functions_mod