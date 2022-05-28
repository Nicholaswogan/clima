
module clima_types
  use clima_const, only: dp, s_str_len
  implicit none
  public
  
  !!!!!!!!!!!!!!!!
  !!! Settings !!!
  !!!!!!!!!!!!!!!!
  
  type :: SettingsOpacity
    character(:), allocatable :: k_method
    integer :: nbins
    
    character(s_str_len), allocatable :: k_distributions(:)
    character(s_str_len), allocatable :: cia(:)
    character(s_str_len), allocatable :: rayleigh(:)
    logical, allocatable :: rayleigh_bool
    character(s_str_len), allocatable :: absorption_xs(:)
    character(s_str_len), allocatable :: photolysis_xs(:)
  end type
  
  type :: ClimaSettings  
    character(:), allocatable :: filename
    
    ! atmosphere-grid
    logical :: atmos_grid_is_present
    real(dp) :: bottom
    real(dp) :: top
    integer :: nz
    
    ! planet
    logical :: planet_is_present
    character(:), allocatable :: back_gas_name
    real(dp) :: P_surf
    real(dp) :: planet_mass
    real(dp) :: planet_radius
    real(dp) :: surface_albedo
    real(dp) :: diurnal_fac
    real(dp) :: solar_zenith
    
    ! optical-properties
    character(s_str_len), allocatable :: species(:)
    type(SettingsOpacity), allocatable :: sol
    type(SettingsOpacity), allocatable :: ir
    
  end type
  
  interface
    module function create_ClimaSettings(filename, err) result(s)
      use fortran_yaml_c, only : parse, error_length
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

end module