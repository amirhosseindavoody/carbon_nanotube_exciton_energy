module cnt_class
	implicit none
	public

	type cnt
		integer :: n_ch,m_ch !chiral vector parameters
		integer :: i_sub !subband index used in exciton energy calculation
	
		!Geometrical properties
		real*8, dimension(2) :: a1,a2,b1,b2,ch_vec,t_vec,aCC_vec
		real*8 :: len_ch,radius
		integer :: Nu !number of graphene unit cells in cnt unit cell.
		integer :: nr !length of cnt in terms of its unit cell.
		real*8, dimension(:,:), allocatable :: posA,posB,posAA,posBB,posAB,posBA
		real*8, dimension(:,:), allocatable :: posA3, posB3
		real*8, dimension(:,:,:), allocatable :: pos2d, pos3d

		!Length and location of cnt for calculating the resonance energy transfer rate
		real*8 :: Length
		real*8 :: center_position

		!Environment properties
		real*8 :: kappa !this is the dielectric factor that accounts for core electron and environment screening
		real*8 :: Ckappa !this is the factor that is multiplied to a kappa_coeff to yield kappa. This is varient under different environments.
		real*8 :: kappa_coeff !this is the scaling factor for calculating kappa and is different for each different cnt chirality.
		
		!Reciprocal lattice properties
		integer :: nkg
		real*8 :: dk
		real*8, dimension(2) :: K1, K2
		
		!CNT band structure properties
		integer, dimension(:), allocatable :: min_sub
		integer :: ikc_max, ikc_min, ik_max, ik_min, iKcm_max, iKcm_min, ik_high, ik_low, ikr_high, ikr_low, iq_max, iq_min
		
		!Dielectric function
		real*8, dimension(:,:), allocatable :: eps_q
		complex*16, dimension(:,:,:,:), allocatable :: v_FT ! v_FT(mu,q,n,m) stores the Fourier transform of the Coulomb potential at the wavevector determined by band index "mu" and wavenumber "q" for atoms of type A (n or m = 1) or type B (n or m = 2)

		!CNT self energy and tight binding coefficients
		real*8, dimension(:,:,:), allocatable :: Ek!Tight-binding energy , Ek(mu,k,n) stores the tight-binding energy of the band "mu" with wavevector "k" in conduction band (n=1) or the valence band (n=2).
		real*8, dimension(:,:,:), allocatable :: Sk!Self-energy 
		complex*16, dimension(:,:,:), allocatable :: Cc,Cv !Tight-binding wavefunction coefficients

		!Exciton wavefunction and energies
		real*8, dimension(:,:), allocatable :: Ex_A1, Ex0_A2, Ex1_A2 !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable :: Psi_A1, Psi0_A2, Psi1_A2 !the first index is ikr, the scond index is the subband, the third index is iKcm

		!Target exciton wavefunction and energies
		real*8, dimension(:,:), allocatable :: Ex_t !the first index is subband, the second index is iKcm
		complex*16, dimension(:,:,:), allocatable :: Psi_t !the first index is ikr, the scond index is the subband, the third index is iKcm
		character(len=20) :: targetExcitonType !this is the type of target exciton which should be one this options: Ex_A1, Ex0_A2, Ex1_A2
		real*8 :: ex_symmetry

		!number of exciton bands below free-electron free-hole energy level
		integer :: nX
		real*8 :: E_th
		real*8 :: Kcm_max

		!directory that the CNT information is stored
		character(len=200) :: directory
	
	end type cnt

end module cnt_class