!*******************************************************************************
! This subroutines opens all the files that are used in the simulation
!*******************************************************************************
subroutine fnOpenFiles()
	use fileHandle
	implicit none
  
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
end
    
!*******************************************************************************
! This subroutines closes all the files that are used in the simulation
!*******************************************************************************
subroutine fnCloseFiles()
	use fileHandle
	implicit none
  
	close(fh2)
	close(fh3)
	close(fh4)
	close(fh5)
	close(fh6)
	close(fh7)
	close(fh8)
	close(fh9)
	close(fh10)
	close(fh11)
	close(fh12)
	close(fh13)
	close(fh14)
	close(fh15)
	close(fh16)
	close(fh17)
	close(fh18)
	close(fh19)
  
	return
end

!*******************************************************************************
! This subroutines opens the log file add new log and closes the file
!*******************************************************************************
subroutine fnLogFile()
	use fileHandle
	implicit none
	
	logical :: exist
	
	inquire(file="log.dat",exist=exist)
	if (exist) then
		open(logFile, file="log.dat", status="old", position="append", action="write")
	else
		open(logFile, file="log.dat", status="new", action="write")
	end if
	write(logFile, *) logInput
	close(logFile)

	return
end subroutine