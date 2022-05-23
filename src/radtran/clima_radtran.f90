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
  
  function create_ClimaRadtranIR(datadir, filename, nz, err) result(rad)
    use clima_types, only: ClimaSettings
    use clima_input, only: create_ClimaSettings
    
    character(*), intent(in) :: datadir
    character(*), intent(in) :: filename
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    type(ClimaSettings) :: s
    
    s = create_ClimaSettings(filename, err)
    if (allocated(err)) return
    
    rad = create_ClimaRadtranIR_(datadir, s%species, s, nz, err)
    if (allocated(err)) return
    
  end function
  
  function create_ClimaRadtranIR_(datadir, species_names, s, nz, err) result(rad)
    use clima_radtran_types_create, only: create_OpticalProperties
    use clima_radtran_types, only: OpticalProperties
    use clima_radtran_types, only: FarUVOpticalProperties, SolarOpticalProperties, IROpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    character(s_str_len), intent(in) :: species_names(:)
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    rad%ng = size(species_names)
    rad%species_names = species_names
    
    rad%ir = create_OpticalProperties(datadir, IROpticalProperties, species_names, s%ir, err)
    if (allocated(err)) return
    
  end function
  
end module
  
  