!*******************************************************************************
! This program calculates the excitonic band structure of single wall carbon nanotubes through simple tight-binding method.
! Amirhossein Davoody
! Last modified: 8/27/2015
!*******************************************************************************

program cnt_exciton_energy
	use A_exciton_energy_mod, only: calculate_A_exciton_dispersion
	use cnt_band_structure_mod, only: cnt_band_structure
	use cnt_geometry_mod, only: cnt_geometry
	use dielectric_fn_mod, only: calculate_dielectric_fn
	use E_exciton_energy_mod, only: calculate_E_exciton_dispersion
	use interpolate_mod, only: interpolate_energy, interpolate_dielectric_fn
	use parse_input_file_mod, only: parse_input_file, finalize_output_directory_name
	use self_energy_mod, only: calculate_self_energy
	use write_log_mod, only: writeLog
	
	implicit none
	
	real :: starttime,endtime !time and date of simulation
	character(len=1000) :: logInput

	call CPU_time(starttime)
	
	call parse_input_file()
	
	call cnt_geometry()
	call cnt_band_structure()

	call calculate_dielectric_fn()
	call interpolate_dielectric_fn()

	call calculate_self_energy()
	call interpolate_energy()

	call calculate_A_exciton_dispersion()
	call calculate_E_exciton_dispersion()

	call finalize_output_directory_name()
  
	call CPU_time(endtime)
	write(logInput,'(A, F0.3, A)') "run time = ",endtime-starttime," seconds."
	call writeLog(trim(logInput))

end program cnt_exciton_energy

