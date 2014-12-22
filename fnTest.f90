!*******************************************************************************
! This subroutines is a test function
!*******************************************************************************
subroutine fnTest()
  
  implicit none
  
  complex*16 , dimension(2) :: arr1,arr2,arr3,arr4
  complex*16 :: num1,num2
  integer :: i, j, fh
  complex*16, parameter :: i1=(0,1.d0)
  
  num1=1+i1
  print *,'abs(num1)=',abs(num1)
  
  num2=0+i1
  print *,'abs(num2)=',abs(num2)
  
  arr1(1)= 1.d0/sqrt(2.d0)
  arr1(2)= 1.d0/sqrt(2.d0)*i1
  
  arr2(1)= 1.d0/sqrt(2.d0)
  arr2(2)= -1.d0/sqrt(2.d0)*i1
  
  arr3=dot_product(arr1,arr1)
  
  arr4=dot_product(arr1,arr2)
  
  pause
  
  return
end