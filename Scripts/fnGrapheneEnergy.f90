subroutine fnGrapheneEnergy(E,Cc,Cv,k)
  use comparams, only: e2p, t0, s0, a1, a2, i1
  !use IMSL_libraries
  implicit none
  
  complex*16 :: f_k
  real*8, dimension(2) :: k
  real*8, dimension(2) :: E
  !real*8, dimension(2) :: tmpE
  !real :: error
  !complex*16, dimension(2,2) :: H, tmpC
  complex*16, dimension(2) :: Cv
  complex*16, dimension(2) :: Cc
  
  f_k=exp(i1*dot_product(k,(a1+a2)/3.d0))+exp(i1*dot_product(k,(a1-2.d0*a2)/3.d0))+exp(i1*dot_product(k,(a2-2.d0*a1)/3.d0))
  !H(1,1)=e2p
  !H(1,2)=t0*f_k
  !H(2,1)=t0*conjg(f_k)
  !H(2,2)=e2p
  !
  !call evchf(H,tmpE,tmpC)
  !error=epihf(2,H,tmpE,tmpC)
  !if (error .ge. 1) stop "*** poor performance in calculating graphene energies ***"
  !
  !if (tmpE(1) .gt. tmpE(2)) then
  !    E(1)=tmpE(1)
  !    Cc=tmpC(:,1)
  !    E(2)=tmpE(2)
  !    Cv=tmpC(:,2)
  !else
  !    E(1)=tmpE(2)
  !    Cc=tmpC(:,2)
  !    E(2)=tmpE(1)
  !    Cv=tmpC(:,1)
  !endif
  !
  !Cc=Cc*conjg(Cc(1))/abs(Cc(1))
  !Cv=Cv*conjg(Cv(1))/abs(CV(1))
  
  E(1)=+t0*abs(f_k)
  E(2)=-t0*abs(f_k)
  
  Cc(1)=+1.d0/sqrt(2.d0)
  Cc(2)=+1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
  Cv(1)=+1.d0/sqrt(2.d0)
  Cv(2)=-1.d0/sqrt(2.d0)*conjg(f_k)/abs(f_k)
  
  return
end