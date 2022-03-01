
module clima_types
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_2d
  implicit none
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  type :: ClimaSettings
    
    
  end type
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Optical Properties and thermodynamics !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type :: Ktable
    integer :: ngauss 
    integer :: ntemp
    integer :: npress
    integer :: nwav
    real(dp), allocatable :: weights(:) ! (ngauss)
    real(dp), allocatable :: temp(:) ! (ntemp)
    real(dp), allocatable :: press(:) ! (nwav)
    type(linear_interp_2d), allocatable :: kappa(:,:) ! (ngauss, nwav)
  end type
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  !!!!!!!!!!!!!!!
  !!! Species !!!
  !!!!!!!!!!!!!!!
  
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
    ! K-distribution
    integer :: ng_k
    type(Ktable), allocatable :: k_solar(:)
    type(Ktable), allocatable :: k_ir(:)
    integer, allocatable :: k_sp_inds(:)
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