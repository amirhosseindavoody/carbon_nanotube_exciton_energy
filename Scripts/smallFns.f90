!**********************************************************************************************************************
! This subroutine calculates the greatest common divisor of the arguments na and nb
!**********************************************************************************************************************
subroutine gcd(ngcd,na,nb)
  integer :: na,nb,ngcd
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
end

!**********************************************************************************************************************
! This subroutine calculates the eigen values and eigen vectors of matrix given the size of it (nl) and gives back eigen vectors in A.
!**********************************************************************************************************************
subroutine eig(nl,matrix,A,W)
	implicit none
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
end subroutine