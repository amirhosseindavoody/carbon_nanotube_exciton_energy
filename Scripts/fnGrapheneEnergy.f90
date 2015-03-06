subroutine fnGrapheneEnergy(E,Cc,Cv,k)
	use comparams, only: t0, a1, a2, i1
	implicit none
	
	complex*16 :: f_k
	real*8, dimension(2) :: k
	real*8, dimension(2) :: E
	complex*16, dimension(2) :: Cv
	complex*16, dimension(2) :: Cc
	
	f_k=exp(i1*dot_product(k,(a1+a2)/3.d0))+exp(i1*dot_product(k,(a1-2.d0*a2)/3.d0))+exp(i1*dot_product(k,(a2-2.d0*a1)/3.d0))
	
	E(1)=+t0*abs(f_k)
	E(2)=-t0*abs(f_k)
	
	Cc(1)=+1.d0/sqrt(2.d0)
	Cc(2)=+1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
	Cv(1)=+1.d0/sqrt(2.d0)
	Cv(2)=-1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
	
	return
end subroutine