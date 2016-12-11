!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parse_input_file_mod
	implicit none
	private

	character(len=1000) :: tmp_output_directory

	public :: parse_input_file, finalize_output_directory_name

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
			write(*,'(A)') "Input format ERROR!"
			write(*,'(A)') "Correct input format is: main.exe inputfile.in"
			call exit()
		end if


		call get_command_argument(1,input_filename)
		open(unit=100,file=input_filename,status="old", action="read", iostat=istat)
		if (istat .ne. 0) then
			write(*,*) 
			write(*,'(A, A)') "Unable to read input file:", input_filename
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
						case ('selected_exciton')
							read(value, *) currcnt%selected_exciton_name
						case ('length[nm]')
							read(value, *) currcnt%length
						case ('center_position[nm]')
							read(value, *) currcnt%center_position
						case default
							write(*,'(A)') "ERROR in 'cnt' input arguments!!!"
							write(*,'(A)') "simulation STOPPED!!!"
							call exit()
						end select
					case ('flg')
						select case (trim(label))
						case ('flg_dielectric')
							read(value, *) flg_dielectric
						case default
							write(*,'(A)') "ERROR in 'flg' input arguments!!!"
							write(*,'(A)') "simulation STOPPED!!!"
							call exit()
						end select
					end select
				end if
			else if (ios .gt. 0) then
				write (*,'(A)') "Error in reading input file!"
				call exit()
			end if
		end do
		close(100)

		! calculate kappa based on input parameters
		if ((currcnt%Ckappa .gt. 0.d0 ) .and. (currcnt%kappa_coeff .gt. 0.d0)) then
			currcnt%kappa = currcnt%Ckappa*currcnt%kappa_coeff
		end if

		! create the name of temporary output directory in which the cnt information is saved
		write(tmp_output_directory,'(A, A, I0, A, I0, A, I0, A, I0, A, I0, A, F0.1, A, F0.1, A, I0, A, F0.1, A)') trim(outdir), "r.exciton_", currcnt%n_ch, "_", currcnt%m_ch, "_nkg_", currcnt%nkg, "_dk_ratio_", currcnt%dk_dkx_ratio, "_nr_", currcnt%nr, "_Eth_", currcnt%E_th/eV, "_Kcm_max_", currcnt%Kcm_max*1.d-9, "_sub_", currcnt%i_sub, "_Ckappa_", currcnt%Ckappa, "/"

		! create the output directory which is the name of the final output directory. the name of the output directory will be changed from tmp_output_directory to cnt%directory when the simulation is finished
		write(currcnt%directory,'(A, A, I0, A, I0, A, I0, A, I0, A, I0, A, F0.1, A, F0.1, A, I0, A, F0.1, A)') trim(outdir), "exciton_", currcnt%n_ch, "_", currcnt%m_ch, "_nkg_", currcnt%nkg, "_dk_ratio_", currcnt%dk_dkx_ratio, "_nr_", currcnt%nr, "_Eth_", currcnt%E_th/eV, "_Kcm_max_", currcnt%Kcm_max*1.d-9, "_sub_", currcnt%i_sub, "_Ckappa_", currcnt%Ckappa, "/"

		! create the output directory and change the working director to that one.
		call create_outdir(tmp_output_directory, input_filename)

		! write simulation settings to the log file
		call writeLog(new_line('A')//"Simulation properties *****************")
		write(logInput,"('n_ch = ',I0)") currcnt%n_ch
		call writeLog(trim(logInput))
		write(logInput,"('m_ch = ',I0)") currcnt%m_ch
		call writeLog(trim(logInput))
		write(logInput,"('nkg = ',I0)") currcnt%nkg
		call writeLog(trim(logInput))
		write(logInput,"('dk/dkx = ',I0)") currcnt%dk_dkx_ratio
		call writeLog(trim(logInput))
		write(logInput,"('nr = ',I0)") currcnt%nr
		call writeLog(trim(logInput))
		write(logInput,"('E_th[eV] = ',F0.1)") currcnt%E_th/eV
		call writeLog(trim(logInput))
		write(logInput,"('Kcm_max[1/nm] = ',F0.1)") currcnt%Kcm_max*1.d-9
		call writeLog(trim(logInput))
		write(logInput,"('i_sub = ',I0)") currcnt%i_sub
		call writeLog(trim(logInput))
		write(logInput,"('Ckappa = ',F0.3)") currcnt%Ckappa
		call writeLog(trim(logInput))
		write(logInput,"('kappa_coeff = ',F0.3)") currcnt%kappa_coeff
		call writeLog(trim(logInput))
		write(logInput,"('kappa = ',F0.3)") currcnt%kappa
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
		write(logInput,'("Simulation started at--> DATE=",I0,"/",I0,"/",I0,"  TIME=",I0,":",I0,":",I0)') date, time
		call writeLog(logInput)

	end subroutine create_outdir


	!***************************************************************************
	! -	this subroutine renames the output directory from a temporary name to a
	!	final name that indicates the simulation has run successfully.
	!***************************************************************************
	subroutine finalize_output_directory_name()
		use comparams, only: currcnt

		logical :: folder_exists
		character(len=1000) :: command
		integer :: istat

		! remove the final output directory if it already exists
		folder_exists = .true.
		inquire(file=trim(currcnt%directory)//'/.', exist=folder_exists)

		if (folder_exists) then
			write(command, "(A, A)") "rm -r ", trim(currcnt%directory)
			call system(trim(command))
		end if

		!rename the temporary output directory to the final output directory
		write(command, '(A, A, A, A)') "mv ", trim(tmp_output_directory), " ", trim(currcnt%directory)
		call system(trim(command))

		!change the working directory to the final output directory
		istat=chdir(trim(currcnt%directory))
		if (istat .ne. 0) then
			write(*,'(A)') "Directory did not changed!!!"
			write(*,'(A)') "Simulation stopped!!!"
			call exit()
		end if

	end subroutine finalize_output_directory_name


end module parse_input_file_mod
