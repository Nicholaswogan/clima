
module clima_types
  use clima_const, only: dp, s_str_len
  use clima_saturationdata, only: SaturationData
  implicit none
  public
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!

  type :: SettingsParticleOpacity
    character(:), allocatable :: name
    character(:), allocatable :: dat
  end type
  
  type :: SettingsOpacity
    integer, allocatable :: new_num_k_bins
  
    character(:), allocatable :: k_method
    integer :: nbins
    
    character(s_str_len), allocatable :: k_distributions(:)
    character(s_str_len), allocatable :: cia(:)
    character(s_str_len), allocatable :: rayleigh(:)
    logical, allocatable :: rayleigh_bool
    character(s_str_len), allocatable :: absorption_xs(:)
    logical, allocatable :: photolysis_bool
    character(s_str_len), allocatable :: photolysis_xs(:)
    character(:), allocatable :: water_continuum
    type(SettingsParticleOpacity), allocatable :: particle_xs(:)
  end type
  
  type :: ClimaSettings  
    character(:), allocatable :: filename
    
    ! atmosphere-grid
    logical :: atmos_grid_is_present
    real(dp), allocatable :: bottom
    real(dp), allocatable :: top
    integer :: nz 
    
    ! planet
    logical :: planet_is_present
    character(:), allocatable :: back_gas_name
    real(dp), allocatable :: P_surf
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp), allocatable :: surface_albedo
    real(dp), allocatable :: diurnal_fac
    real(dp), allocatable :: solar_zenith
    
    ! optical-properties
    character(:), allocatable :: wavelength_bins_file
    character(s_str_len), allocatable :: gases(:)
    character(s_str_len), allocatable :: particles(:)
    type(SettingsOpacity), allocatable :: sol
    type(SettingsOpacity), allocatable :: ir
    
  end type
  
  interface
    module function create_ClimaSettings(filename, err) result(s)
      character(*), intent(in) :: filename
      character(:), allocatable, intent(out) :: err
      type(ClimaSettings) :: s
    end function
  end interface
  interface ClimaSettings
    module procedure :: create_ClimaSettings
  end interface
  
  
  type :: AtmosphereFile
    character(:), allocatable :: filename
    integer :: nz
    integer :: nlabels
    character(s_str_len), allocatable :: labels(:)
    real(dp), allocatable :: columns(:,:) ! (size(labels),nz)
  end type
  
  interface
    module function create_AtmosphereFile(atm_file, err) result(atm)
      character(*), intent(in) :: atm_file
      character(:), allocatable, intent(out) :: err
      type(AtmosphereFile) :: atm
    end function
  end interface
  interface AtmosphereFile
    module procedure create_AtmosphereFile
  end interface
  
  ! other
  interface
    module subroutine unpack_atmospherefile(atm, species_names, z, mix, T, P, err)
      type(AtmosphereFile), intent(in) :: atm
      character(*), intent(in) :: species_names(:)
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: mix(:,:)
      real(dp), intent(out) :: T(:)
      real(dp), intent(out) :: P(:)
      character(:), allocatable, intent(out) :: err
    end subroutine
  end interface
  
  !!!!!!!!!!!!!!!
  !!! Species !!!
  !!!!!!!!!!!!!!!
  
  enum, bind(c)
  enumerator :: &
    ShomatePolynomial = 1, &
    Nasa9Polynomial = 2
  end enum
  
  type :: ThermodynamicData
    integer :: dtype
    integer :: ntemps
    real(dp), allocatable :: temps(:)
    real(dp), allocatable :: data(:,:)
  end type
  
  type :: Gas
    character(:), allocatable :: name
    integer, allocatable :: composition(:) ! (natoms)
    real(dp) :: mass
    
    ! thermodynamics
    type(ThermodynamicData) :: thermo
    type(SaturationData), allocatable :: sat
    
  end type
  
  type :: Species
    integer :: natoms
    character(s_str_len), allocatable :: atoms_names(:)
    real(dp), allocatable :: atoms_mass(:)
    
    integer :: ng
    type(Gas), allocatable :: g(:) ! (ng)
  end type
  
  interface
    module function create_Species(filename, err) result(sp)
      character(*), intent(in) :: filename
      character(:), allocatable, intent(out) :: err
      type(Species) :: sp
    end function
  end interface
  interface Species
    module procedure create_Species
  end interface

end module