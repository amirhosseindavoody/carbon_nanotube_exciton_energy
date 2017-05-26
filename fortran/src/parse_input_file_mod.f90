!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Declaration of input parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module parse_input_file_mod
	implicit none
	private

	character(len=1000) :: tmp_output_directory

	public :: parse_input_file

contains
	!**************************************************************************************************************************
	! parse the input file
	!**************************************************************************************************************************

	subroutine parse_input_file()
		use comparams, only: currcnt
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
		currcnt%name='cnt1'
		currcnt%n=10
		currcnt%m=0
		currcnt%nr=100
		currcnt%ex%i_sub=1
		
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
						case ('name')
							currcnt%name = trim(value)
						case ('n')
							read(value, *) currcnt%n
						case ('m')
							read(value, *) currcnt%m
						case ('length[unit cells]')
							read(value, *) currcnt%nr
						case ('exciton_subband')
							read(value, *) currcnt%ex%i_sub
						case ('selected_exciton')
							read(value, *) currcnt%ex%name
						case default
							write(*,'(A)') "ERROR in 'cnt' input arguments!!!"
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

		! create the output directory which is the name of the final output directory. the name of the output directory will be changed from tmp_output_directory to cnt%directory when the simulation is finished
		! write(currcnt%directory,'(A, A, I0, A, I0, A, I0, A, I0, A, I0, A, F0.1, A, F0.1, A, I0, A, F0.1, A)') trim(outdir), "exciton_", currcnt%n_ch, "_", currcnt%m_ch, "_nkg_", currcnt%nkg, "_dk_ratio_", currcnt%dk_dkx_ratio, "_nr_", currcnt%nr, "_Eth_", currcnt%E_th/eV, "_Kcm_max_", currcnt%Kcm_max*1.d-9, "_sub_", currcnt%i_sub, "_Ckappa_", currcnt%Ckappa, "/"
		currcnt%directory = outdir



		! create the output directory and change the working director to that one.
		call create_outdir(outdir, input_filename)

		! write simulation settings to the log file
		call writeLog(new_line('A')//"Simulation properties *****************")
		write(logInput,"('n = ',I0)") currcnt%n
		call writeLog(trim(logInput))
		write(logInput,"('m = ',I0)") currcnt%m
		call writeLog(trim(logInput))
		write(logInput,"('nr = ',I0)") currcnt%nr
		call writeLog(trim(logInput))
		write(logInput,"('exciton_subband = ',I0)") currcnt%ex%i_sub
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
		write(command,'("rm -rf ",A)') trim(outdir) !remove the directory if it already exists
		call system(trim(command))
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


	! !***************************************************************************
	! ! -	this subroutine renames the output directory from a temporary name to a
	! !	final name that indicates the simulation has run successfully.
	! !***************************************************************************
	! subroutine finalize_output_directory_name()
	! 	use comparams, only: currcnt

	! 	logical :: folder_exists
	! 	character(len=1000) :: command
	! 	integer :: istat

	! 	! remove the final output directory if it already exists
	! 	folder_exists = .true.
	! 	inquire(file=trim(currcnt%directory)//'/.', exist=folder_exists)

	! 	if (folder_exists) then
	! 		write(command, "(A, A)") "rm -r ", trim(currcnt%directory)
	! 		call system(trim(command))
	! 	end if

	! 	!rename the temporary output directory to the final output directory
	! 	write(command, '(A, A, A, A)') "mv ", trim(tmp_output_directory), " ", trim(currcnt%directory)
	! 	call system(trim(command))

	! 	!change the working directory to the final output directory
	! 	istat=chdir(trim(currcnt%directory))
	! 	if (istat .ne. 0) then
	! 		write(*,'(A)') "Directory did not changed!!!"
	! 		write(*,'(A)') "Simulation stopped!!!"
	! 		call exit()
	! 	end if

	! end subroutine finalize_output_directory_name

end module parse_input_file_mod
