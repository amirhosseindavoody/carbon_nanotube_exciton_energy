!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program CNT_Exciton
  use comparams
  implicit none
  
  !call fnTest
  
  call fnInput
  
  call fnOpenFiles
  
  call fnSaveSimInfo
  
  call fnPhysConst
  
  call fnGeomProp
  
  call fnCNTband

  call fnDielectric
  
  call fnSelfEnergy
  
  call fnExcitonDispersion
  
  call fnSaveMisc
  
  call fnCloseFiles
  

  print *,'Finish!!!!'
  
endprogram

