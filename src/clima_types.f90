
module clima_types
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
  implicit none
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  type :: ClimaSettings
    
  end type

  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Optical Properties  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: Ktable
    integer :: sp_ind
    integer :: ngauss 
    integer :: ntemp
    integer :: npress
    integer :: nwav
    real(dp), allocatable :: weights(:) ! (ngauss)
    real(dp), allocatable :: temp(:) ! (ntemp)
    real(dp), allocatable :: press(:) ! (nwav)
    type(linear_interp_2d), allocatable :: kappa(:,:) ! (ngauss, nwav)
  end type
  
  type :: CIAtable
    integer :: sp_inds(2)
    integer :: ntemp
    integer :: nwav
    real(dp), allocatable :: temp(:) ! (ntemp)
    type(linear_interp_1d), allocatable :: kappa(:) ! (nwav) ! [cm^-1*(cm^3/molecules)^-2]
  end type
  
  type :: XsectionData
    integer :: sp_ind
    integer :: ntemp
    real(dp), allocatable :: temp(:) ! (ntemp)
    real(dp), allocatable :: xs(:,:) ! (ntemp, nw)
  end type
  
  !!!!!!!!!!!!!!!
  !!! Species !!!
  !!!!!!!!!!!!!!!
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  type :: Species
    character(:), allocatable :: name
    integer, allocatable :: composition(:) ! (natoms)
    real(dp) :: mass
    
    ! thermodynamics
    type(ThermodynamicData) :: thermo
    
  end type

  type :: ClimaData
    
    integer :: natoms
    character(s_str_len), allocatable :: atoms_names(:)
    real(dp), allocatable :: atoms_mass(:)
    
    integer :: ng
    type(Species), allocatable :: sp(:)
    
    !!! Optical properties !!!
    ! K-distributions (e.g. H2O)
    integer :: nk_sol
    type(Ktable), allocatable :: k_sol(:)
    integer :: nk_ir
    type(Ktable), allocatable :: k_ir(:)
    ! T-dependent CIA coefficients (e.g. H2-H2)
    integer :: ncia_sol
    type(CIAtable), allocatable :: cia_sol(:)
    integer :: ncia_ir
    type(CIAtable), allocatable :: cia_ir(:)
    ! Cross sections (e.g. O3 photolysis)
    integer :: nxs_sol
    type(XsectionData), allocatable :: xs_sol(:)
    integer :: nxs_ir
    type(XsectionData), allocatable :: xs_ir(:)
    ! Rayleigh Scattering
    integer :: nray
    real(dp), allocatable :: sigray(:,:) ! (nray, nw)
    integer, allocatable :: ray_sp_inds(:) ! species number of rayleigh species
    
  end type
  
  type :: ClimaVars
    
    character(:), allocatable :: data_dir
    integer :: nz
    real(dp), allocatable :: mix(:)
  
  end type
  
  type :: ClimaWrk
    ! work variables
  end type
  
end module