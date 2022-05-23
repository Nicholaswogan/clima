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
  
contains
  
  
  function create_ClimaRadtranIR(datadir, s, nz, err) result(rad)
    use clima_const, only: c_light
    use clima_radtran_types, only: OpticalProperties
    use clima_types, only: ClimaSettings
    
    character(*), intent(in) :: datadir
    type(ClimaSettings), intent(in) :: s
    integer, intent(in) :: nz
    character(:), allocatable, intent(out) :: err
    
    type(ClimaRadtranIR) :: rad
    
    
  end function
  
end module
  
  