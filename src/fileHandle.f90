!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are the unit numbers for the input/output files in the program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fileHandle
	implicit none
	
	character(len=200) :: logInput
	integer, parameter :: logFile = 10
	integer :: fh1,fh2,fh3,fh4,fh5,fh6,fh7,fh8,fh9,fh10,fh11,fh12,fh13,fh14,fh15,fh16,fh17,fh18,fh19

contains
	
	!*******************************************************************************
	! This subroutines opens all the files that are used in the simulation
	!*******************************************************************************
	
	subroutine fnOpenFiles()
	  
		fh2=11
		open(unit=fh2,file='posA.dat',status="unknown")
	  
		fh3=12
		open(unit=fh3,file='posB.dat',status="unknown")
	  
		fh4=13
		open(unit=fh4,file='CondBand.dat',status="unknown")
	  
		fh5=14
		open(unit=fh5,file='ValeBand.dat',status="unknown")
	  
		fh6=15
		open(unit=fh6,file='CondBand_Sub.dat',status="unknown")
	  
		fh7=16
		open(unit=fh7,file='ValeBand_Sub.dat',status="unknown")
	  
		fh8=17
		open(unit=fh8,file='PI_q.dat',status="unknown")
	  
		fh9=18
		open(unit=fh9,file='v_FT.dat',status="unknown")
	  
		fh10=19
		open(unit=fh10,file='eps_q.dat',status="unknown")
	  
		fh11=20
		open(unit=fh11,file='CondSelfEnergy_Sub.dat',status="unknown")
	  
		fh12=21
		open(unit=fh12,file='ValeSelfEnergy_Sub.dat',status="unknown")
	  
		fh13=22
		open(unit=fh13,file='Ex_A1.dat',status="unknown")
	  
		fh14=23
		open(unit=fh14,file='Ex0_A2.dat',status="unknown")
	  
		fh15=24
		open(unit=fh15,file='Ex1_A2.dat',status="unknown")
	  
		fh16=25
		open(unit=fh16,file='Psi_A1.dat',status="unknown")
	  
		fh17=26
		open(unit=fh17,file='Psi0_A2.dat',status="unknown")
	  
		fh18=27
		open(unit=fh18,file='Psi1_A2.dat',status="unknown")
	  
		fh19=28
		open(unit=fh19,file='miscellaneous.dat',status="unknown")

		return
	end subroutine fnOpenFiles

end module fileHandle