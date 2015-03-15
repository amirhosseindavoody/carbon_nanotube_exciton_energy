!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 3/17/2014
!*******************************************************************************

program CNT_Exciton
	use comparams
	use fileHandle
	implicit none
	
	call CPU_time(starttime)
	
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
  
	call CPU_time(endtime)
	write(logInput,'("Run time = ",f10.3," seconds.")'),endtime-starttime
	call fnLogFile()
  
endprogram

