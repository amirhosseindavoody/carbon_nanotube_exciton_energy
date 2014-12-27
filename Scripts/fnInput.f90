!*******************************************************************************
! This subroutines interprets the input arguments of the simulation
!*******************************************************************************
subroutine fnInput()
  use comparams
  implicit none
  
  integer i,count,flg_tmp
  character*20 :: buffer
  
  ! set the simulation variables to default values
  n_ch=10
  m_ch=0
  nkg=501
  nr=200
  E_th=1.5d0
  Kcm_max=1.5d9
  flg_dielectric=.true.  !when .true. dielectric function is calculated, when .false. dielectric function is read from file.
  i_sub=1
  
  ! update the simulation variables according to input variables
  count=nargs()
  i=1
  
  do while (i .le. count-1)
    call getarg(i,buffer)
    if (buffer .eq. 'ch') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) n_ch
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) m_ch
    elseif (buffer .eq. 'nkg') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) nkg
    elseif (buffer .eq. 'nr') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) nr
    elseif (buffer .eq. 'E_th') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) E_th      
    elseif (buffer .eq. 'Kcm_max') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) Kcm_max
      Kcm_max=Kcm_max*1.d9
    elseif (buffer .eq. 'flg_dielectric') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) flg_tmp
      if (flg_tmp .eq. 1) then
        flg_dielectric=.true.
      elseif (flg_tmp .eq. 0) then
        flg_dielectric=.false.
      else
        print *, "ERROR in input argument flg_dielectric!"
        pause
        stop
      endif
    elseif (buffer .eq. 'i_sub') then
      i=i+1
      call getarg(i,buffer)
      read(buffer,*) i_sub 
    else
        print *, "ERROR in input arguments!"
        pause
        stop
    end if
    i=i+1
  end do

  return
end