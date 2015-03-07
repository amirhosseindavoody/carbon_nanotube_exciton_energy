!*******************************************************************************
! This subroutines interprets the input arguments of the simulation
!*******************************************************************************
subroutine fnInput()
  use comparams
  use fileHandle
  implicit none
  
  integer, dimension(3) :: date, time
  integer i,count,flg_tmp
  character(len=20) :: buffer
  
  ! get time and date of start of simulation
  call idate(date)
  call itime(time)
  
  ! set the simulation variables to default values
  n_ch=10
  m_ch=0
  nkg=0501
  nr=200
  E_th=1.5d0
  Kcm_max=1.5d9
  flg_dielectric=.true.  !when .true. dielectric function is calculated, when .false. dielectric function is read from file.
  i_sub=1
  kappa=2.d0
  
  ! update the simulation variables according to input variables
  count=command_argument_count()
  i=1
  
  do while (i .le. count-1)
    call get_command_argument(i,buffer)
    if (buffer .eq. 'ch') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) n_ch
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) m_ch
    elseif (buffer .eq. 'nkg') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) nkg
    elseif (buffer .eq. 'nr') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) nr
    elseif (buffer .eq. 'E_th') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) E_th      
    elseif (buffer .eq. 'Kcm_max') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) Kcm_max
      Kcm_max=Kcm_max*1.d9
    elseif (buffer .eq. 'flg_dielectric') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) flg_tmp
      if (flg_tmp .eq. 1) then
        flg_dielectric=.true.
      elseif (flg_tmp .eq. 0) then
        flg_dielectric=.false.
      else
        write(logInput,*) "ERROR in input argument flg_dielectric!"
		call fnLogFile()
        read (*,*)
        stop
      endif
    elseif (buffer .eq. 'i_sub') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) i_sub
	elseif (buffer .eq. 'kappa') then
      i=i+1
      call get_command_argument(i,buffer)
      read(buffer,*) kappa
    else
        write(logInput,*) "ERROR in input arguments!"
		call fnLogFile()
        read(*,*)
        stop
    end if
    i=i+1
  end do
  
  ! write simulation inputs to the log file
  write(logInput,*) "Simulation started at: DATE=", date,"    TIME=", time
  call fnLogFile()
  write(logInput,*) "SIMULATION PARAMETERS"
  call fnLogFile()
  write(logInput,*) "n_ch=",n_ch
  call fnLogFile()
  write(logInput,*) "m_ch=",m_ch
  call fnLogFile()
  write(logInput,*) "nkg=",nkg
  call fnLogFile()
  write(logInput,*) "nr=",nr
  call fnLogFile()
  write(logInput,*) "E_th=",E_th
  call fnLogFile()
  write(logInput,*) "Kcm_max=",Kcm_max
  call fnLogFile()
  write(logInput,*) "flg_dielectric=",flg_dielectric
  call fnLogFile()
  write(logInput,*) "i_sub=",i_sub
  call fnLogFile()
  write(logInput,*) "kappa=",kappa
  call fnLogFile()

  return
end