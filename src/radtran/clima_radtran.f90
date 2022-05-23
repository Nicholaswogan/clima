module clima_radtran
  use clima_const, only: dp, s_str_len
  use clima_radtran_types, only: OpticalProperties, RadiateXSWrk, RadiateZWrk
  implicit none
  
  ! initialize with
  ! settings.yaml, sun.txt, and data dir, species names
  ! 
  ! inputs are
  ! surface_albedo, u0, T, densities, dz
  ! ouptuts are
  ! fluxes (fup, fdn)
  
  type :: ClimaRadtranIR
  
    integer :: ng
    character(s_str_len), allocatable :: species_names(:) ! (ng) copy of species names
  
    !!! Optical properties !!!
    type(OpticalProperties) :: ir
    
    type(RadiateXSWrk) :: rx_ir
    type(RadiateZWrk) :: rz_ir
  
  end type
  
  type, extends(ClimaRadtranIR) :: ClimaRadtran

    !!! Optical properties !!!
    type(OpticalProperties) :: sol
  
    real(dp) :: diurnal_fac = 0.5_dp
    real(dp), allocatable :: photons_sol(:) ! (nw) mW/m2/Hz in each bin  
    
    type(RadiateXSWrk) :: rx_sol
    type(RadiateZWrk) :: rz_sol
    
  end type
  
  interface ClimaRadtranIR
    module procedure :: create_ClimaRadtranIR
  end interface
  
contains
  
  function create_ClimaRadtranIR(datadir, settings_f, nz, err) result(rad)
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: settings_f
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    type(ClimaSettings) :: s
    
    s = ClimaSettings(settings_f, err)
    if (allocated(err)) return
    
    if (.not. allocated(s%species)) then
      err = '"'//settings_f//'/optical-properties/species" does not exist'
      return
    endif
    
    rad = create_ClimaRadtranIR_(datadir, s%species, s, nz, err)
    if (allocated(err)) return
    
  end function
  
  function create_ClimaRadtranIR_(datadir, species_names, s, nz, err) result(rad)
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(s_str_len), intent(in) :: species_names(:)
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    if (nz < 1) then
      err = '"nz" can not be less than 1.'
      return
    endif
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    
    if (.not. allocated(s%ir)) then
      err = '"'//s%filename//'/optical-properties/ir" does not exist.'
      return
    endif
    rad%ir = OpticalProperties(datadir, IROpticalProperties, species_names, s%ir, err)
    if (allocated(err)) return
    
    rad%rx_ir = RadiateXSWrk(rad%ir, nz)
    rad%rz_ir = RadiateZWrk(nz)

  end function

  
end module
  
  