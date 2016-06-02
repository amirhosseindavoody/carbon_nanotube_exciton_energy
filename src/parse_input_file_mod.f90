!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parse_input_file_mod
	implicit none
	private
	public :: parse_input_file

contains
	!**************************************************************************************************************************
	! parse the input file
	!**************************************************************************************************************************

	subroutine parse_input_file()
		use comparams, only: currcnt, flg_dielectric
		use physical_constant_mod, only: eV
		use write_log_mod, only: writeLog

		character(len=1000) :: input_filename
		character(len=1000) :: buffer, command, label, value
		character(len=1000) :: outdir
		character(len=1000) :: logInput
		integer :: istat=0
		integer :: ios=0
		integer :: pos_comma=0, pos_equal=0

		if (command_argument_count() .ne. 1) then
			write(*,*) "Input format ERROR!"
			write(*,*) "Correct input format is: main.exe inputfile.in"
			call exit()
		end if


		call get_command_argument(1,input_filename)
		open(unit=100,file=input_filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,*) ""
			write(*,*) "Unable to read input file:", input_filename
			call exit()
		end if

		! set the simulation variables to default values
		currcnt%n_ch=10
		currcnt%m_ch=0
		currcnt%nkg=0501
		currcnt%nr=200
		currcnt%E_th=1.5d0
		currcnt%Kcm_max=1.5d9
		currcnt%i_sub=1
		currcnt%kappa=2.d0
		currcnt%Ckappa=0.d0
		currcnt%kappa_coeff=0.d0
		flg_dielectric=.true.  !when .true. dielectric function is calculated, when .false. dielectric function is read from file.

		do while (ios == 0)
			read (100,'(A)',iostat=ios) buffer
			if (ios == 0) then
				if (buffer(1:1) .ne. '#') then
					pos_comma = scan(buffer,',')
					pos_equal = scan(buffer,'=')
					command = adjustl(buffer(1:pos_comma-1))
					label = adjustl(buffer(pos_comma+1:pos_equal-1))
					value = adjustl(buffer(pos_equal+1:))

					! set the target cnt
					select case (trim(command))
					case('directory')
						select case (trim(label))
						case ('output')
							outdir = trim(value)
						case default
							write(*,*) "ERROR in 'directory' input arguments!!!"
							write(*,*) "simulation STOPPED!!!"
							call exit()
						end select
					case ('cnt')
						select case (trim(label))
						case ('n_ch')
							read(value, *) currcnt%n_ch
						case ('m_ch')
							read(value, *) currcnt%m_ch
						case ('nkg')
							read(value, *) currcnt%nkg
						case ('dk/dkx')
							read(value, *) currcnt%dk_dkx_ratio
						case ('nr')
							read(value, *) currcnt%nr
						case ('E_th[eV]')
							read(value, *) currcnt%E_th
							currcnt%E_th = currcnt%E_th * eV
						case ('Kcm_max[1/nm]')
							read(value, *) currcnt%Kcm_max
							currcnt%Kcm_max = currcnt%Kcm_max * 1.d9
						case ('i_sub')
							read(value, *) currcnt%i_sub
						case ('Ckappa')
							read(value, *) currcnt%Ckappa
						case ('kappa_coeff')
							read(value, *) currcnt%kappa_coeff
						case ('target_exciton_type')
							read(value, *) currcnt%targetExcitonType
						case ('length[nm]')
							read(value, *) currcnt%length
						case ('center_position[nm]')
							read(value, *) currcnt%center_position
						case default
							write(*,*) "ERROR in 'cnt' input arguments!!!"
							write(*,*) "simulation STOPPED!!!"
							call exit()
						end select
					case ('flg')
						select case (trim(label))
						case ('flg_dielectric')
							read(value, *) flg_dielectric
						case default
							write(*,*) "ERROR in 'flg' input arguments!!!"
							write(*,*) "simulation STOPPED!!!"
							call exit()
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,*) "Error in reading input file!"
				call exit()
			end if
		end do
		close(100)

		! calculate kappa based on input parameters
		if ((currcnt%Ckappa .gt. 0.d0 ) .and. (currcnt%kappa_coeff .gt. 0.d0)) then
			currcnt%kappa = currcnt%Ckappa*currcnt%kappa_coeff
		end if

		! create the output directory in which the cnt information is saved
		write(currcnt%directory,"( A, 'CNT(', I2.2, ',', I2.2, ')-nkg(', I4.4, ')-nr(', I4.4, ')-E_th(', F3.1, ')-Kcm_max(', F3.1, ')-i_sub(', I1.1, ')-Ckappa(', F3.1, ')/' )") trim(outdir), currcnt%n_ch, currcnt%m_ch, currcnt%nkg, currcnt%nr, currcnt%E_th/eV, currcnt%Kcm_max*1.d-9, currcnt%i_sub, currcnt%Ckappa

		! create the output directory and change the working director to that one.
		call create_outdir(currcnt%directory, input_filename)

		! write simulation settings to the log file
		call writeLog(new_line('A')//"Simulation properties *****************")
		write(logInput,"('n_ch = ',I2.2)") currcnt%n_ch
		call writeLog(trim(logInput))
		write(logInput,"('m_ch = ',I2.2)") currcnt%m_ch
		call writeLog(trim(logInput))
		write(logInput,"('nkg = ',I4.4)") currcnt%nkg
		call writeLog(trim(logInput))
		write(logInput,"('dk/dkx = ',I4.4)") currcnt%dk_dkx_ratio
		call writeLog(trim(logInput))
		write(logInput,"('nr = ',I4.4)") currcnt%nr
		call writeLog(trim(logInput))
		write(logInput,"('E_th[eV] = ',F3.1)") currcnt%E_th/eV
		call writeLog(trim(logInput))
		write(logInput,"('Kcm_max[1/nm] = ',F3.1)") currcnt%Kcm_max*1.d-9
		call writeLog(trim(logInput))
		write(logInput,"('i_sub = ',I1.1)") currcnt%i_sub
		call writeLog(trim(logInput))
		write(logInput,"('Ckappa = ',F3.1)") currcnt%Ckappa
		call writeLog(trim(logInput))
		write(logInput,"('kappa_coeff = ',F5.3)") currcnt%kappa_coeff
		call writeLog(trim(logInput))
		write(logInput,"('kappa = ',F5.3)") currcnt%kappa
		call writeLog(trim(logInput))

	end subroutine parse_input_file

	!*******************************************************************************
	! This subroutines changes the working directory to the output directory for saving files
	!*******************************************************************************

	subroutine create_outdir(outdir, input_filename)
		use write_log_mod, only: writeLog

		character(len=*), intent(in) :: outdir
		character(len=*), intent(in) :: input_filename
		integer :: istat=0
		character(len=1000) :: command
		integer, dimension(3) :: date, time
		character(len=1000) :: logInput

		! specifiy the output directory
		write(command,'("rm -rf ''",A,"''")') trim(outdir) !remove the directory if it already exists
		call system(command)
		write(command,'("mkdir ''",A,"''")') trim(outdir) !create the directory again
		call system(trim(command))

		! copy the input files to the output directory
		write(command,'(A)') "cp '"//trim(input_filename)//"' '"//trim(outdir)//"'"
		call system(trim(command))

		istat=chdir(trim(outdir))
		if (istat .ne. 0) then
			write(*,'(A)') "Directory did not changed:"
			write(*,'(A)') trim(outdir)
			write(*,'(A)') "Simulation stopped!!!"
			call exit()
		end if

		! get time and date of start of simulation
		call idate(date)
		call itime(time)

		! write simulation inputs to the log file
		write(logInput,'("Simulation started at--> DATE=",I2.2,"/",I2.2,"/",I4.4,"  TIME=",I2.2,":",I2.2,":",I2.2)') date, time
		call writeLog(logInput)

		return
	end subroutine create_outdir


end module parse_input_file_mod
