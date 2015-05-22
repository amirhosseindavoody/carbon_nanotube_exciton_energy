module write_log_mod
	implicit none
	private
	public :: writeLog

contains
	!*******************************************************************************
	! This subroutines opens the log file add new log and closes the file
	!*******************************************************************************

	subroutine writeLog(message)
		character(len=*) :: message
		logical :: flgexist
		integer :: logFile = 10

		inquire(file="log.dat",exist=flgexist)
		if (flgexist) then
			open(logFile, file="log.dat", status="old", position="append", action="write")
		else
			open(logFile, file="log.dat", status="new", action="write")
		end if
		write(logFile, *) message
		close(logFile)

		return
	end subroutine writeLog
	
end module write_log_mod