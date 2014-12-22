!**********************************************************************************************************************
! This subroutines interprets the input arguments of the simulation
!**********************************************************************************************************************
subroutine fnGeomProp()
  use comparams
  implicit none
  integer :: tmpi,i,j,k
  real*8 :: p1,p2,p,q,tmpr
  real*8 :: ndc,mdc,Nur
  real*8 :: cosTh, sinTh
  real*8, dimension(2,2) :: Rot
  real*8, dimension(2) :: tmp_vec
  ! unit vectors and reciprocal lattice vectors************************************************************************
  ndc=n_ch
  mdc=m_ch
  a1=(/dsqrt(3.d0)/2.d0*a_l, 1.d0/2.d0*a_l/)
  a2=(/dsqrt(3.d0)/2.d0*a_l, -1.d0/2.d0*a_l/)
  b1=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, +1.d0*2.d0*pi/a_l/)
  b2=(/1.d0/dsqrt(3.d0)*2.d0*pi/a_l, -1.d0*2.d0*pi/a_l/)
  
  aCC_vec=1.d0/3.d0*(a1+a2)
  
  ! calculate chirality and translational vectors of CNT unit cell.****************************************************
  ch_vec=dble(n_ch)*a1+dble(m_ch)*a2

  len_ch=a_l*dsqrt(dble(n_ch)**2+dble(m_ch)**2+dble(n_ch)*dble(m_ch))
  radius=len_ch/2.d0/pi
  
  call gcd(dR,2*n_ch+m_ch,2*m_ch+n_ch)
  
  t1=(2.d0*dble(m_ch)+dble(n_ch))/dble(dR)
  t2=-(2.d0*dble(n_ch)+dble(m_ch))/dble(dR)
  
  t_vec=t1*a1+t2*a2
  
  Nu=2*(n_ch**2+m_ch**2+n_ch*m_ch)/dR

  p1=max((dble(Nu)/dble(n_ch)+1.d0/t1)/(dble(m_ch)/dble(n_ch)-t2/t1),(1.d0/dble(n_ch)+1.d0/t1)/(dble(m_ch)/dble(n_ch)-t2/t1))
  p2=min((dble(Nu)/dble(n_ch)+1.d0/t1)/(dble(m_ch)/dble(n_ch)-t2/t1),(1.d0/dble(n_ch)+1.d0/t1)/(dble(m_ch)/dble(n_ch)-t2/t1))
  
  do i=ceiling(p2),floor(p1)
      p=i
      q=dble(t2)/dble(t1)*p+1.d0/dble(t1)
      if (q .eq. ceiling(q)) then
        exit
      else if (i==floor(p1)) then
        write(fh1,10) "MC not found"
        stop
      end if
  end do
  
  call gcd(tmpi,int(abs(p)),int(abs(q)))
  p=p/dble(tmpi)
  q=q/dble(tmpi)
  MC=dble(m_ch)*p-dble(n_ch)*q
  
  ! rotate basis vectors so that ch_vec is along x-axis
  cosTh=ch_vec(1)/norm2(ch_vec)
  sinTh=ch_vec(2)/norm2(ch_vec)
  Rot=reshape((/ cosTh, -sinTh , sinTh, cosTh /), (/2,2/))
  ch_vec=matmul(Rot,ch_vec)
  t_vec=matmul(Rot,t_vec)
  a1=matmul(Rot,a1)
  a2=matmul(Rot,a2)
  b1=matmul(Rot,b1)
  b2=matmul(Rot,b2)
  aCC_vec=matmul(Rot,aCC_vec)
  
  ! calculate reciprocal lattice of CNT.*******************************************************************************
  dk=norm2(b1)/(dble(nkg)-1.d0)
  K1=(-t2*b1+t1*b2)/(dble(Nu))
  K2=(dble(m_ch)*b1-dble(n_ch)*b2)/dble(Nu)
  K2=K2/norm2(K2)
  
  ! calculate coordinates of atoms in the unwarped CNT unit cell.******************************************************
  allocate(posA(Nu,2))
  allocate(posB(Nu,2))
    
  k=0
  do i=0,t1+n_ch
    do j=t2,m_ch
      if ((dble(t2)/dble(t1)*i .le. j) .and. (dble(m_ch)/dble(n_ch)*i .ge. j) .and. (dble(t2)/dble(t1)*(i-n_ch) .gt. (j-m_ch)) .and. (dble(m_ch)/dble(n_ch)*(i-t1) .lt. (j-t2))) then
        k=k+1
        posA(k,1)=dble(i)*a1(1)+dble(j)*a2(1)
        posA(k,2)=dble(i)*a1(2)+dble(j)*a2(2)
        posB(k,1)=posA(k,1)+aCC_vec(1)
        posB(k,2)=posA(k,2)+aCC_vec(2)
        
        if (posA(k,1) .gt. ch_vec(1)) posA(k,1)=posA(k,1)-ch_vec(1);
        if (posA(k,1) .lt. 0) posA(k,1)=posA(k,1)+ch_vec(1);
        if (posA(k,2) .gt. t_vec(2)) posA(k,2)=posA(k,2)-t_vec(2);
        if (posA(k,2) .lt. 0) posA(k,2)=posA(k,2)+t_vec(2);
        
        if (posB(k,1) .gt. ch_vec(1)) posB(k,1)=posB(k,1)-ch_vec(1);
        if (posB(k,1) .lt. 0) posB(k,1)=posB(k,1)+ch_vec(1);
        if (posB(k,2) .gt. t_vec(2)) posB(k,2)=posB(k,2)-t_vec(2);
        if (posB(k,2) .lt. 0) posB(k,2)=posB(k,2)+t_vec(2);
        
        write(fh2,14) posA(k,1), posA(k,2)
        write(fh3,14) posB(k,1), posB(k,2)
        
      endif
    enddo
  enddo
    
  if (k .ne. Nu) stop "*** Error in calculating atom positions ***"
  
  ! calculate distances between atoms in a warped CNT unit cell.*******************************************************
  allocate(posAA(Nu,2))
  allocate(posAB(Nu,2))
  allocate(posBA(Nu,2))
  allocate(posBB(Nu,2))
  
  do i=1,Nu
    posAA(i,:)=posA(i,:)-posA(1,:)
    posAB(i,:)=posA(i,:)-posB(1,:)
    posBA(i,:)=posB(i,:)-posA(1,:)
    posBB(i,:)=posB(i,:)-posB(1,:)
    if (posAA(i,1) .gt. ch_vec(1)/2.d0) posAA(i,1)=posAA(i,1)-ch_vec(1)
    if (posAB(i,1) .gt. ch_vec(1)/2.d0) posAB(i,1)=posAB(i,1)-ch_vec(1)
    if (posBA(i,1) .gt. ch_vec(1)/2.d0) posBA(i,1)=posBA(i,1)-ch_vec(1)
    if (posBB(i,1) .gt. ch_vec(1)/2.d0) posBB(i,1)=posBB(i,1)-ch_vec(1)
    
    !write(fh2,14) posAA(i,1), posAA(i,2)
    !write(fh3,14) posBA(i,1), posBA(i,2)
    
  end do
  
  ! write down important informations into the output file.************************************************************
  write(fh1,10) "GEOMETRICAL PROPERTIES"
  write(fh1,11) "a1=",a1(1), a1(2)
  write(fh1,11) "a2=",a2(1), a2(2)
  write(fh1,11) "b1=",b1(1), b1(2)
  write(fh1,11) "b2=",b2(1), b2(2)
  write(fh1,11) "aCC_vec=",aCC_vec(1), aCC_vec(2)
  write(fh1,11) "ch_vec=",ch_vec(1), ch_vec(2)
  write(fh1,11) "t_vec=",t_vec(1), t_vec(2)
  write(fh1,12) "len_ch=",len_ch
  write(fh1,13) "Nu=",Nu
  write(fh1,13) "MC=",MC
  

10 FORMAT (A100)
11 FORMAT (A10,E16.8,E16.8)   
12 FORMAT (A10,E16.8)
13 FORMAT (A10,I5)
14 FORMAT (E16.8,E16.8)  
   
  return
end