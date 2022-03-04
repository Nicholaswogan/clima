
module clima_types
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
  implicit none
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  
  type :: SettingsOpacity
    character(:), allocatable :: name
    character(:), allocatable :: type
    character(s_str_len), allocatable :: species(:)
    character(:), allocatable :: solar_filename
    character(:), allocatable :: ir_filename
  end type
  
  type :: ClimaSettings
    
    type(SettingsOpacity), allocatable :: op(:)
    
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
    real(dp), allocatable :: log10P(:) ! (nwav)
    type(linear_interp_2d), allocatable :: log10kappa(:,:) ! (ngauss, nwav)
  end type
  
  type :: CIAtable
    integer :: sp_inds(2)
    integer :: ntemp
    integer :: nwav
    real(dp), allocatable :: temp(:) ! (ntemp)
    type(linear_interp_1d), allocatable :: log10kappa(:) ! (nwav) ! [cm^-1*(cm^3/molecules)^-2]

  end type
  
  type :: Xsection
    integer :: sp_ind
    integer :: ntemp
    real(dp), allocatable :: temp(:) ! (ntemp)
    ! if ntemp == 1, use xs
    real(dp), allocatable :: xs(:) ! (nw) 
    ! else, use xs_i
    type(linear_interp_1d), allocatable :: xs_i(:) ! (nw) 
  end type
  
  enum, bind(c)
    enumerator :: SolarOpticalProperties, IROpticalProperties
  end enum
  
  type :: OpticalProperties
    integer :: type
    integer :: nw
    real(dp), allocatable :: wavenums(:)
    
    ! K-distributions (e.g. H2O)
    integer :: nk
    type(Ktable), allocatable :: k(:)
    ! T-dependent CIA coefficients (e.g. H2-H2)
    integer :: ncia
    type(CIAtable), allocatable :: cia(:)
    ! Cross sections (e.g. O3 photolysis)
    integer :: nxs
    type(Xsection), allocatable :: xs(:)
    ! Rayleigh Scattering
    integer :: nray
    real(dp), allocatable :: sigray(:,:) ! (nray, nw)
    integer, allocatable :: ray_sp_inds(:) ! species number of rayleigh species
    
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
    character(s_str_len), allocatable :: species_names(:) ! copy of species name
    type(Species), allocatable :: sp(:)
    
    !!! Optical properties !!!
    type(OpticalProperties) :: sol
    type(OpticalProperties) :: ir
    
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