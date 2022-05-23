
module clima_types
  use iso_c_binding
  use clima_const, only: dp, s_str_len
  use linear_interpolation_module, only: linear_interp_1d, linear_interp_2d
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
    type(SettingsOpacity) :: sol
    type(SettingsOpacity) :: ir
    
  end type
  
  type :: AtmosphereFile
    integer :: nz
    integer :: nlabels
    character(s_str_len), allocatable :: labels(:)
    real(dp), allocatable :: columns(:,:) ! (size(labels),nz)
  end type

end module